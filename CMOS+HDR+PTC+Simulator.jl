
using Distributions, FixedPointNumbers, Gadfly, Images, ImageView, Colors, StatsBase, DataFrames, IterTools

# type definitions to shorten function parameter lists
type AnalogChain
  em_gain::Float64      # EM gain when present.  Set to 0 to bypass
  cv_gain::Float64      # conversion gain in uV/e-
  read_noise::Float64   # read/source follower noise in input e- equivalent
  pga_gain::Float64     # PGA gain in uV/uV
  non_linear::Float64   # non linearity as a % deviation at FWC
  adc_offset::Float64   # offset into the ADC in uV
  bit_depth::Float64    # ADC bit depth
  min_volts::Float64    # ADC minimum input voltage
  max_volts::Float64    # ADC maximum input voltage
end

type ImageDimensions
  x_dim::Int            # columns
  y_dim::Int            # rows
  z_dim::Int            # frames
end

# returns an image stack representing detected photo-electrons
function photons2electrons(photon_number, QE, id::ImageDimensions)
  rand(Poisson(QE*photon_number), id.x_dim, id.y_dim, id.z_dim)
end

# em register transfer function
function em_register(photoelectrons, ac::AnalogChain)
  try
    x_dim, y_dim, z_dim = size(photoelectrons)
  end

  read_noise = rand(Poisson(signal), x_dim, y_dim, z_dim)
  # offset(uV) + cvg(uV/e) (signal(e-) + read_noise(e-))
  ac.adc_offset + ac.pga_gain*(ac.cv_gain*(photoelectrons + read_noise))
end

# analog transfer functions
function analog_chain_v(photoelectrons, ac::AnalogChain)
  try
    x_dim, y_dim, z_dim = size(photoelectrons)
  end

  read_noise = rand(Normal(0, ac.read_noise), x_dim, y_dim, z_dim)
  # offset(uV) + cvg(uV/e) (signal(e-) + read_noise(e-))
  ac.adc_offset + ac.pga_gain*(ac.cv_gain*(photoelectrons + read_noise))
end

# (residual) fixed pattern noise model: assumes a random distribution between columns and rows
# every pixel gets a random amount of gain/offset variation.
function fpn(z, col_offset, col_gain, row_offset, row_gain, idim::ImageDimensions)

  # co = column offset, cg = column gain
  covec = rand(Normal(0, col_offset), idim.y_dim)
  cgvec = rand(Normal(0, col_gain), idim.y_dim)
  rovec = reshape(rand(Normal(0, row_offset), idim.x_dim), 1, idim.x_dim)
  rgvec = reshape(rand(Normal(0, row_gain), idim.x_dim), 1, idim.x_dim)

  z = covec .+ z
  z = cgvec .* z
  z = rovec .+ z
  z = rgvec .* z
end

# convert PGA voltages to digital numbers
function analog_to_digital_v(values, ac::AnalogChain)
  cvalues = round.(values*(2^ac.bit_depth-1)/(ac.max_volts - ac.min_volts))
  clamp.(cvalues, 0, 2^ac.bit_depth-1)
end

# high signal non-linear gain transform
function nonlinear(s_uv, nl2_term)
   nl2_term * s_uv^2
end

# dual output HDR Gain Combiner logic
# gc_p: Noise weighted averaging
function gc_p(g1, g2, g1_min_signal, g1_max_signal)
  if g1 < g1_min_signal
    out = g1
  elseif g1 <= g1_max_signal
    out = 2 * g1 * g2 / (g1 + g2 )
  elseif g1 >= g1_max_signal
    out = g2
  else
    out = 0
  end
  out
end

# gc_a: Straight averaging gain combiner
function gc_a(g1, g2, g1_min_signal, g1_max_signal)
  if g1 < g1_min_signal
    out = g1
  elseif g1 <= g1_max_signal
    out = (g1 + g2 )/2
  elseif g1 >= g1_max_signal
    out = g2
  else
    out = 0
  end
  out
end

# dual gain HDR
function gain_combiner_e_v(g1_img_e, g2_img_e, gcf, g1_min::Int, g1_max::Int)
  gcf.(g1_img_e, g2_img_e, g1_min, g1_max)
end

function gain_combiner_dn_v(g1_img_e, g2_img_e, g2_fullwell_e::Int, gcf, out_bits::Int, out_offset::Int, g1_min::Int, g1_max::Int)
  out_offset .+ round.(((2^out_bits-1)/g2_fullwell_e) .* gcf.(g1_img_e, g2_img_e, g1_min, g1_max))
end


function ptc1(min_s, max_s, step, ac::AnalogChain, idim::ImageDimensions )
  p = []
  r = []
  r2 = []

  for i in imap(x -> x^2, range(min_s, step, Int(round((max_s-min_s)/step))))
    z = analog_chain_v(i, ac, idim)
    q = analog_to_digital_v(z, ac)

    append!( p, mean(mean(q,3))) #perhaps most probable/frequent
    append!( r, mean(std(q,3)))
    append!( r2, mean(var(q,3)))
  end
  p = p - round(ac.adc_offset*(2^ac.bit_depth-1)/(ac.max_volts - ac.min_volts))
  q = DataFrame(Signal_ADU=p, Noise_ADU=r, Variance_ADU2=r2)
  q
end

# generate a dual gain HDR PTC plot, return 3 DataFrames
function ptc2(min_s, max_s, step, ac1::AnalogChain, ac2::AnalogChain, id::ImageDimensions)
  a1m = []
  a1s = []
  a1v = []
  a2m = []
  a2s = []
  a2v = []
  cm = []
  cs = []
  cv = []
  const out_bits = 16
  const out_offset_dn = 102
  const cmb_range_min = 1700
  const cmb_range_max = 1900

  # recast gain and offset to a more convienient form
  num_steps = Int(round((max_s-min_s)/step))
  s4o = (ac1.adc_offset*(2^ac1.bit_depth-1)/(ac1.max_volts - ac1.min_volts))
  s4g = (ac1.max_volts-ac1.min_volts)/(ac1.cv_gain * ac1.pga_gain * (2^ac1.bit_depth - 1))
  q4o = (ac2.adc_offset*(2^ac2.bit_depth-1)/(ac2.max_volts - ac2.min_volts))
  q4g = (ac2.max_volts-ac2.min_volts)/(ac2.cv_gain * ac2.pga_gain * (2^ac2.bit_depth - 1))

  for i in imap(x -> x, range(min_s, step, num_steps))
    s0 = photons2electrons(i, 0.8, id)
    s1 = s0 #em_multiplier(s0,ac1.em_gain)
    s2 = analog_chain_v(s1, ac1)
    s3 += nonlinear.(s2, 0)
    s3 = analog_to_digital_v(s3, ac1)
    s4 = clamp.((s4 .- s4o),0, 65535) .* s4g

    # For fun - include a bit of non linearity in the low gain channel
    q0 = s0  #(same input)
    q1 = em_multiplier(q0, ac2.em_gain)
    q2 = analog_chain_v(q1, ac2)
    q3 += nonlinear.(q2, -0.05/400000)
    q4 = analog_to_digital_v(q3, ac2)
    q5 = clamp.((q4 .- q4o), 0, 65536) .* q4g

    q3 = gain_combiner_e_v(q1, q2, gc_p, cmb_range_min, cmb_range_max)

    # stats on each image stack
    append!( a1m, mean(mean(q1,3))) # perhaps most probable/frequent
    append!( a1s, mean(std(q1,3)))
    append!( a1v, mean(var(q1,3)))

    append!( a2m, mean(mean(q2,3))) # perhaps most probable/frequent
    append!( a2s, mean(std(q2,3)))
    append!( a2v, mean(var(q2,3)))

    append!( cm, mean(mean(q3,3))) # perhaps most probable/frequent
    append!( cs, mean(std(q3,3)))
    append!( cv, mean(var(q3,3)))
  end

  # PTC DataFrames for inputs plus combined output
  a1_df = DataFrame(Signal_E=a1m, Noise_E=a1s, Variance_E2=a1v)
  a2_df = DataFrame(Signal_E=a2m, Noise_E=a2s, Variance_E2=a2v)
  c_df = DataFrame(Signal_E=cm, Noise_E=cs, Variance_E2=cv)
  a1_df, a2_df, c_df
end


# example with gain combiner
# (em gain, conversion gain (uV/e), read noise(electrons), PGA, non linearity, offset(uV), adc bits, min volts, max volts)
  a1 = AnalogChain( 1, 40, 1.5, 3, 0, 20000, 11, 0, 400000)
  a2 = AnalogChain( 1, 10, 15, 2, 0, 20000, 11, 0, 400000)
  id12 = ImageDimensions( 100, 100, 100)

# covert offset in uV to DN, system gains in e/DN, FW in electrons
  of1_DN = (a1.adc_offset*(2^a1.bit_depth-1)/(a1.max_volts - a1.min_volts))
  of2_DN = (a2.adc_offset*(2^a2.bit_depth-1)/(a2.max_volts - a2.min_volts))

  sg1 = 1/(a1.cv_gain*a1.pga_gain*(2^a1.bit_depth-1)/(a1.max_volts - a1.min_volts))
  sg2 = 1/(a2.cv_gain*a2.pga_gain*(2^a2.bit_depth-1)/(a2.max_volts - a2.min_volts))

  fw1_electrons = Int(round((2^a1.bit_depth-1) * sg1))
  fw2_electrons = Int(round((2^a2.bit_depth-1) * sg2))

  sgc = fw2_electrons/(2^16-1)

(fw1_electrons, fw2_electrons)

#(min signal(e), max signal(e), step(e), analog chain 1, analog chain 2, image dimensions )
  ptc_a1, ptc_a2, ptc_c = ptc2(0, 21000, 100, a1, a2, id12)

#pretty plots
rick_morty = Theme(point_size = 3pt, point_size_max = 3pt)

Gadfly.push_theme(rick_morty)

p4 = plot(layer(ptc_c, x="Signal_E", y="Noise_E", Geom.point, order=1, Theme(default_color=colorant"Red")),
     layer(ptc_a1, x="Signal_E", y="Noise_E", Geom.line, order=2, Theme(default_color=colorant"Green")),
     layer(ptc_a2, x="Signal_E", y="Noise_E", Geom.line, order=3, Theme(default_color=colorant"Blue")),
     Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
     Guide.xlabel("Mean Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve: sCMOS HDR"))

p5 = plot(layer(ptc_c, x="Signal_E", y="Variance_E2", Geom.point, order=1, Theme(default_color=colorant"Red")),
     layer(ptc_a1, x="Signal_E", y="Variance_E2", Geom.line, order=2, Theme(default_color=colorant"Green")),
     layer(ptc_a2, x="Signal_E", y="Variance_E2", Geom.line, order=3, Theme(default_color=colorant"Blue")),
     Guide.xlabel("Mean Signal(e)"), Guide.ylabel("Variance(e2)"), Guide.title("Mean Variance Curve: sCMOS HDR"))

#display image of a 100e- signal - averaged across 100 frames
a1 = AnalogChain( 40, 1.5, 3, 20000, 11, 0, 400000)
id12 = ImageDimensions( 100, 100, 100)
z1 = analog_chain_v(100, a1, id12)
q1 = analog_to_digital_v(z1, a1)
colorview(Gray, mean(q1,3)[:,:,1]/512)

em_img = rand(Gamma(10, 100-1+1/10),100,100)
colorview(Gray, em_img/2000)
mean(em_img)
q = DataFrame(ADU = vec(em_img))
plot(q, x="ADU", Geom.histogram)
