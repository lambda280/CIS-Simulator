# PMQ CMOS Camera Simulator
# Version 0.1
#

# run one time only: add these packages to the standard distribution
  Pkg.add("Distributions")
  Pkg.add("FixedPointNumbers")
  Pkg.add("Gadfly")
  Pkg.add("ImageMagick")
  Pkg.add("Images")
  Pkg.add("ImageView")
  Pkg.add("Colors")
  Pkg.add("StatsBase")
  Pkg.add("DataFrames")
  Pkg.add("IterTools")

# One time: add the Jupyter/IJulia for web browser based sessions
Pkg.add("IJulia")

# To launch an IJulia session - run the following from here or a julia terminal prompt
using IJulia; notebook(detached=true)

# run the following each time Julia is started - imports packages into the current session
using Distributions, FixedPointNumbers, Gadfly, Images, ImageView, Colors, StatsBase, DataFrames, IterTools

# type definitions to shorten function parameter lists
type AnalogChain
  cv_gain::Float64      # conversion gain in uV/e-
  read_noise::Float64   # read/source follower noise in input e- equivalent
  pga_gain::Float64     # PGA gain in uV/uV
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

# analog chain: assumes a Poisson distribution of signal, Gaussian for read noise.
#   to do: dark current is not simulated
#   to do: distribute read noise across the various components
function analog_chain_v(photoelectron_number, ac::AnalogChain, id::ImageDimensions)
  ac.adc_offset + ac.pga_gain*(ac.cv_gain*(rand(Poisson(photoelectron_number), id.x_dim, id.y_dim, id.z_dim) + rand(Normal(0, ac.read_noise), id.x_dim, id.y_dim, id.z_dim)))
end

# (residual) fixed pattern noise: assumes a random distribution
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

# dual output HDR logic
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

#generate a single channel PTC plot, return a DataFrame of (mean-bias, std_dev) values
function ptc(min_s, max_s, step, ac::AnalogChain, idim::ImageDimensions )
  p = []
  r = []
  r2 = []

  for i in imap(x -> x^0.5, range(min_s, step, Int(round((max_s-min_s)/step))))
    z = analog_chain_v(i, ac, idim)
    q = analog_to_digital(z, ac)

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
  const cmb_range_min = 1300
  const cmb_range_max = 1500

  num_steps = Int(round((max_s-min_s)/step))
  q1o = (ac1.adc_offset*(2^ac1.bit_depth-1)/(ac1.max_volts - ac1.min_volts))
  q1g = (ac1.max_volts-ac1.min_volts)/(ac1.cv_gain * ac1.pga_gain * (2^ac1.bit_depth - 1))
  q2o = (ac2.adc_offset*(2^ac2.bit_depth-1)/(ac2.max_volts - ac2.min_volts))
  q2g = (ac2.max_volts-ac2.min_volts)/(ac2.cv_gain * ac2.pga_gain * (2^ac2.bit_depth - 1))

  for i in imap(x -> sqrt(x), range(min_s, step, num_steps))
    z1 = analog_chain_v(i, ac1, id)
    q1 = analog_to_digital_v(z1, ac1)
    q1 = clamp.((q1 .- q1o),0, (2^ac1.bit_depth-1)) .* q1g

    z2 = analog_chain_v(i, ac2, id)
    #z2 += nonlinear.(z2, -0.05/400000)
    q2 = analog_to_digital_v(z2, ac2)
    q2 = clamp.((q2 .- q2o), 0, (2^ac2.bit_depth-1)) .* q2g

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
# (conversion gain, read noise(electrons), PGA, offset(uV), adc bits, min volts, max volts)
  a1 = AnalogChain( 40, 1.5, 3, 20000, 11, 0, 400000)
  a2 = AnalogChain( 10, 20, 2, 20000, 11, 0, 400000)
  id12 = ImageDimensions( 100, 100, 100)

# offset in electrons and DN, system gains, FW in electrons
  of1_DN = (a1.adc_offset*(2^a1.bit_depth-1)/(a1.max_volts - a1.min_volts))
  of2_DN = (a2.adc_offset*(2^a2.bit_depth-1)/(a2.max_volts - a2.min_volts))

  sg1 = 1/(a1.cv_gain*a1.pga_gain*(2^a1.bit_depth-1)/(a1.max_volts - a1.min_volts))
  sg2 = 1/(a2.cv_gain*a2.pga_gain*(2^a2.bit_depth-1)/(a2.max_volts - a2.min_volts))

  fw1_electrons = Int(round((2^a1.bit_depth-1) * sg1))
  fw2_electrons = Int(round((2^a2.bit_depth-1) * sg2))

  sgc = fw2_electrons/(2^16-1)

#(min signal(e), max signal(e), step(e), analog chain 1, analog chain 2, image dimensions )
  ptc_a1, ptc_a2, ptc_c = ptc2(0, 21000, 100, a1, a2, id12)

p1 = plot(ptc_c, x="Signal_E", y="Noise_E", Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
            Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve CG"))

p2 = plot(ptc_a1, x="Signal_E", y="Noise_E", Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
            Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve HG"))

x_marks = [0,1,2,3]
y_marks = [0,1,2,3]
p3 = plot(ptc_a2, x="Signal_E", y="Noise_E", Scale.x_log10(minvalue=1), Scale.y_log10(minvalue=10),
            Guide.xticks(ticks=x_marks),
            Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve LG"))

p2a = plot(ptc_a1, x="Signal_E", y="Noise_E", Coord.cartesian(xmin=0, ymin=0), Geom.line, Geom.point,
          Scale.x_log10, Scale.y_log10,
          Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve"))

p2b = plot(ptc_a2, x="Signal_E", y="Variance_E2", Coord.cartesian, Geom.line, Geom.point,
          Scale.x(minvalue=10), Scale.y(minvalue=10),
          Guide.xlabel("Signal(e)"), Guide.ylabel("Variance(e^2)"), Guide.title("Mean Variance Curve"))

p4 = plot(layer(ptc_c, x="Signal_E", y="Noise_E", Geom.point, order=1),
     layer(ptc_a1, x="Signal_E", y="Noise_E", Geom.line, order=2),
     layer(ptc_a2, x="Signal_E", y="Noise_E", Geom.line, order=3),
     Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
     Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve"))

p5 = plot(layer(ptc_c, x="Signal_E", y="Noise_E", Geom.point, order=1),
        layer(ptc_a1, x="Signal_E", y="Noise_E", Geom.line, order=2),
        layer(ptc_a2, x="Signal_E", y="Noise_E", Geom.line, order=3),
        Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
        Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve"))

# statistics through a Z stack -> array
  mean(q,3)
  std(q,3)

# full stack statistics -> value
  mean(mean(q,3))
  mean(std(q,3))

# display an image
  colorview(Gray, mean(z,3)[:,:,1]/512)

# plot histogram
  q = DataFrame(ADU = vec(z))
  plot(q, x="ADU", Geom.histogram)

# example of FPN
  A = ones(1000,1000,10)
  imaged = ImageDimensions(1000,1000,10)
  B = fpn(A, 10, 1.05, 0.1, 0.01, imaged)
  colorview(Gray, (B[:,:,2]-minimum(B))/(maximum(B)-minimum(B)))
