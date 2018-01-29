#~~~~ Copyright Roper Technologies DBA Photometrics/QImaging ~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~ Camera Simulator V0.0 1/28/2018 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using Distributions, GLM, FixedPointNumbers, Gadfly, Images, ImageView, Colors, StatsBase, DataFrames, IterTools

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~ type definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   - immutable types are read-only after instantiation
#   - needed: support for independent saturation in pixel, SF, ADC, etc.
#   - needed: additional noise params for each stage aka SF, ADC, etc.

immutable ImageDimensions
  x_dim::Int32            # columns
  y_dim::Int32            # rows
  z_dim::Int32            # frames
end

immutable PixelProperties
  qe::Float32             # Quantum Efficiency
  x_dim::Float32          # pixel X (um)
  y_dim::Float32          # pixel Y (um)
  i_dark::Float32         # dark current in pA/um2
end

immutable AnalogChain
  em_gain::Float64        # EM gain when present.  Set to 0 to bypass
  em_stages::Int64        # number of stages in EM register
  cv_gain::Float64        # sf conversion gain in uV/e-
  read_noise::Float64     # read noise in input referred e- equivalent
  pga_gain::Float64       # PGA gain in uV/uV
  non_linearity::Float64  # non linearity as a % deviation at FWC
  offset::Float64         # offset in uV used to create th ADC bias offset
end

immutable FPN
  col_offset::Float32
  col_gain::Float32
  row_offset::Float32
  row_gain::Float32
end

immutable ADC
  bit_depth::Int64        # ADC bit depth
  min_volts::Float64      # ADC minimum input voltage
  max_volts::Float64      # ADC maximum input voltage
  bad_range               # A range type which ADC malfunctions
  excess_noise::Float64   # excess noise as percentage of signal
end

immutable ADCCombiner
  bits::Int64        # = 16
  offset_dn::Int64   # = 100
  min::Int64         # = 1000
  max::Int64         # = 1100
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~ optical transfer functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~  electronic transfer functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Issue: The one-pixel approach to drawing random numbers is slow.
#   Using rand(distribution, array) is faster as long as distribution exists in stats module
#   This is not the case for EMCCD excess noise, RTN, or even Idark
#
#   Solutions:
#   - the stats package could be extended to include the needed distributions
#   - where possible create pool of random numbers to draw from inside transfer functions

# returns Poisson distributed random number representing detected photo-electrons
function photon_noise(pixel_flux::Float64, pixel::PixelProperties)
  rand(Poisson(pixel.qe*pixel_flux))
end

# this is a stub
function dark_noise(temperature::Float64, pixel::PixelProperties)
  rand(Poisson(pixel.i_dark))
end

# EMCCD EM register transfer function
#   - A functional form for the distribution function would be faster
function em_register(input_e::Int64, ac::AnalogChain)
  if ac.em_gain <= 0 return input_e end
  accumulator = input_e
  for i in range(1,ac.em_stages)
    accumulator += rand(Binomial(accumulator,ac.em_gain))
  end
  accumulator
end

# returns a Normal distributed random number representing read noise
  # read noise is an input referred amalgamation of several noise factors
  # - move to voltage domain and redestribute across the approprate
  #   stages of the analog chain
function read_noise(ac::AnalogChain)
  rand(Normal(0, ac.read_noise))
end

# stub: returns <distribution> distributed random number representing rtn noise
  # can add noise power spectrum and timing dependencies
function rtn_noise(ac::AnalogChain)
  rand(Normal(0, ac.read_noise))
end

# analog chain transfer functions
  # can add further noise models with modeling CDS process
function analog_chain_v(photoelectrons::Int64, ac::AnalogChain)
  # offset(uV) + cvg(uV/e-) [signal(e-) + read_noise(e-)]
  ac.offset + ac.pga_gain*(ac.cv_gain*n)+ac.pga_gain*ac.cv_gain*read_noise(ac)
end

# fixed pattern noise
  # - the row/col nature breaks the "by pixel" design paradigm
function fpn(z, fp::FPN, i::ImageDimensions)
  # nomenclature: co = column offset, cg = column gain
  covec = fp.col_offset > 0 ? rand(Normal(0, fp.col_offset), i.y_dim) : 0
  cgvec = fp.col_gain > 0 ? rand(Normal(0, fp.col_gain), i.y_dim) : 1
  rovec = fp.row_offset > 0 ? reshape(rand(Normal(0, fp.row_offset), i.x_dim), 1, i.x_dim) : 0
  rgvec = fp.row_gain > 0 ? reshape(rand(Normal(0, fp.row_gain), i.x_dim), 1, i.x_dim) : 1

  z = covec .+ z #< apply column offset
  z = cgvec .* z #< apply column gain
  z = rovec .+ z #< apply row offset
  z = rgvec .* z #< apply row gain
end

# convert PGA output to digital numbers
  # ideally use fixed point math via "FixedPointNumbers"
function analog_to_digital_v(value, adc::ADC)
  cvalue = value*(2^adc.bit_depth-1)/(adc.max_volts - adc.min_volts)
  cvalue = round(((cvalue >= adc.bad_range[1]) && (cvalue <= adc.bad_range[end])) ? rand(Normal(cvalue,adc.excess_noise)) : cvalue)
  clamp(cvalue, 0, 2^adc.bit_depth-1)
end

# signal response nonlinearity transfer function
function nonlinear(signal, ac::AnalogChain)
   signal + ac.non_linearity * signal^2
end

# dual output HDR Gain Combiner logic functions
#   ideally convert to use fixed point math

# gc_p: Noise weighted averaging based on Poisson assumption
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

# dual output HDR Gain Combiner logic
# gc_p_approx: Weighted averaging based on constants in cross over region (variance)
function gc_w(g1, g2, w1, w2, g1_min_signal, g1_max_signal)
  if g1 < g1_min_signal
    out = g1
  elseif g1 <= g1_max_signal
    out = g1/w1 + g2/w2
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
#   this stub is for clarity - can be removed
function gain_combiner_e_e(g1_e::Int32, g2_e::Int32, gcf, g1_min::Int32, g1_max::Int32)
  gcf(g1_e, g2_e, g1_min, g1_max)
end

# same as above, but rescaled into DN and digital offset restored
function gain_combiner_e_dn(g1_e, g2_e, gcf, g2_fullwell_e::Int, out_bits::Int, out_offset::Int, g1_min::Int, g1_max::Int)
  out_offset + Int(round(((2^out_bits-1)/g2_fullwell_e) * gcf(g1_e, g2_e, g1_min, g1_max)))
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~ image analysis functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create a simulated image
function image_sim(signal, pixel::PixelProperties, ac::AnalogChain, fp_n::FPN, adc::ADC, id::ImageDimensions, pattern=ones)
  img = pattern(Float64,id.x_dim, id.y_dim, id.z_dim)

  s0_p = signal .* img
  s1_e = photon_noise.(s0_p, pixel)
  s1_em_e = em_register.(s1_e, ac)
  s2_v = analog_chain_v.(s1_em_e, ac)
  s3_v = fpn(s2_v, fp_n, id)
  s4_v = nonlinear.(s3_v, ac)             #< should fpn be part of analog_chain_v?
  s5_dn = analog_to_digital_v.(s4_v, adc)
end

# Autogeneration of "exposure levels" for pretty PTC plots
function ptc_exposures(min_signal::Int64, max_signal::Int64, steps::Int64)
  step_size = log10(max_signal + 1)/steps
  exposures = imap(x -> (10^x)-1, range(min_signal, step_size, steps))
end


# Photon Transfer and Mean Variance
#   - optional/default params must be at end of argument list
function ptc(pixel::PixelProperties, ac::AnalogChain, fp_noise::FPN, adc::ADC, id::ImageDimensions, pattern_f=ones, expose_f=ptc_exposures)
  a0m_e = Float64[]
  a0s_e = Float64[]
  a0m_dn = Float64[]
  a0v_dn = Float64[]

  a1m_e = Float64[]
  a1s_e = Float64[]
  a1m_dn = Float64[]
  a1v_dn = Float64[]

  # recast offset and gain to a more convienient form
  offset = Int(round(ac.offset*(2^adc.bit_depth-1)/(adc.max_volts - adc.min_volts)))
  gain = (adc.max_volts-adc.min_volts)/(ac.cv_gain * ac.pga_gain * (2^adc.bit_depth - 1))
  fwc = Int(round((2^adc.bit_depth-1) * gain))

  # "autogeneration" of exposure times/signal levels for PTC
  if expose_f == ptc_exposures
    exposures = expose_f(0, fwc, 100)
  #elseif <some sort of test if a function or list was passed in expose_f>
  #   try/except wrapper
  else
    exposures = [0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000]
  end

  # noiseless image "pattern" to use for PTC/MV.  Default is the "ones" function to simulate a flat field
  img = pattern_f(id.x_dim, id.y_dim, id.z_dim)

  em_net_gain = ac.em_gain <= 0 ? 1 : mean(em_register.(img, ac))

  # simulate image sequence taken at multiple exposure times/signal levels
  #  - default signal levels, "exposures", are equally spaced on log graph (better for PTC than MV)
  for i in exposures
    s0_p = i .* img
    s1_e = photon_noise.(s0_p, pixel)
    s1_em_e = em_register.(s1_e, ac)
    s2_v = analog_chain_v.(s1_em_e, ac)
    s3_v = fpn(s2_v, fp_noise, id)
    s4_v = nonlinear.(s3_v, ac)
    s5_dn = analog_to_digital_v.(s4_v, adc)

    #convert DN back to electrons and remove offset for the PTC vector/graph
    s5_dn_no = (s5_dn .- offset) ./ em_net_gain
    s5_e_no =  ((s5_dn .- offset) .* gain) ./ em_net_gain

    # calculate image statistics
    append!( a0m_e, mean(mean(s1_e,3))) # perhaps most probable/frequent
    append!( a0s_e, mean(std(s1_e,3)))

    append!( a0m_dn, mean(mean(s5_dn_no,3)))
    append!( a0v_dn, mean(var(s5_dn_no,3)))

    append!( a1m_e, mean(mean(s5_e_no,3))) # perhaps most probable/frequent
    append!( a1s_e, mean(std(s5_e_no,3)))

    append!( a1m_dn, mean(mean(s5_dn_no,3)))
    append!( a1v_dn, mean(var(s5_dn_no,3)))
  end

  # output image statistics as tuple of DataFrames for use in mv/ptc plots
  s1_df = DataFrame(Signal_E=a0m_e, Noise_E=a0s_e, Signal_DN=a0m_dn, Variance_DN=a0v_dn)
  s5_df = DataFrame(Signal_E=a1m_e, Noise_E=a1s_e, Signal_DN=a1m_dn, Variance_DN=a1v_dn)
  s1_df, s5_df
end

# generate a dual gain HDR PTC plot, return 3 DataFrames
function ptc(pattern=ones, expose_f=ptc_exposures, ac1::AnalogChain, ac2::AnalogChain, adc::ADC, combiner::ADCCombiner, id::ImageDimensions, pixel::PixelProperties)
  a1m = []
  a1s = []
  a1v = []
  a2m = []
  a2s = []
  a2v = []
  cm = []
  cs = []
  cv = []

  # recast offset and gain to a more convienient form
  s4o = Int(round(ac1.adc_offset*(2^ac1.bit_depth-1)/(ac1.max_adc_volts - ac1.min_adc_volts)))
  s4g = (ac1.max_adc_volts-ac1.min_adc_volts)/(ac1.cv_gain * ac1.pga_gain * (2^ac1.bit_depth - 1))
  q4o = Int(round(ac2.adc_offset*(2^ac2.bit_depth-1)/(ac2.max_adc_volts - ac2.min_adc_volts)))
  q4g = (ac2.max_adc_volts-ac2.min_adc_volts)/(ac2.cv_gain * ac2.pga_gain * (2^ac2.bit_depth - 1))
  fwc = Int(round((2^ac2.bit_depth-1) * q4g))
  sg = fwc/(2^16-1)
  step_size = log10(fwc + 1)/steps

  # "generate" is a function that defines object patterns
  #    - default is a flat field
  img = generate(Int16,id.x_dim, id.y_dim, id.z_dim)

  # simulate images at multiple exposure times/signal levels
  #  - signal levels are equally spaced on log graph (better for PTC than MV)
  for i in imap(x -> 10^x-1, range(min_s, step_size, steps))
    s0 = i .* img
    s1 = photon_noise.(s0, pixel)
    #s1 = em_register.(s1, ac1)
    s2 = analog_chain_v.(s1, ac1)
    s2 += nonlinear.(s2, ac1)
    s4 = analog_to_digital_v.(s2, ac1)
    s5 = (s4 .- s4o) .* s4g

    q1 = s1 #< both channels see the same signal in e-, different read noise however
    q2 = analog_chain_v.(q1, ac2)
    q2 += nonlinear.(q2, ac2)
    q4 = analog_to_digital_v.(q2, ac2)
    q5 = (q4 .- q4o) .* q4g

    z3 = gain_combiner_e_dn.(s5, q5, gc_a, fwc, Int(16), Int(0), cmb_range_min, cmb_range_max)
    # convert back to electrons
    z5 = sg .* z3

    # stats on each image stack
    append!( a1m, mean(mean(s5,3))) # perhaps most probable/frequent
    append!( a1s, mean(std(s5,3)))
    append!( a1v, mean(var(s5,3)))

    append!( a2m, mean(mean(q5,3))) # perhaps most probable/frequent
    append!( a2s, mean(std(q5,3)))
    append!( a2v, mean(var(q5,3)))

    append!( cm, mean(mean(z5,3))) # perhaps most probable/frequent
    append!( cs, mean(std(z5,3)))
    append!( cv, mean(var(z5,3)))
  end

  # output image statistics as tuple of DataFrames for use in mv/ptc plots
  a1_df = DataFrame(Signal_E=a1m, Noise_E=a1s, Variance_E2=a1v)
  a2_df = DataFrame(Signal_E=a2m, Noise_E=a2s, Variance_E2=a2v)
  c_df =  DataFrame(Signal_E=cm,  Noise_E=cs,  Variance_E2=cv)
  a1_df, a2_df, c_df
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~ Convienience Functions  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DEFAULT_SAVE_WIN = "C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\"
DEFAULT_SAVE_MAC = "/Users/rob/Dropbox/Simulator/Graphs/"

function ptc_plot(ptc_vecs, title="PTC", fit=false, save=true)
  ticks_x = [0:1:6 ;]
  ticks_y = [0:0.2:3 ;]
  layers = []

  for i in ptc_vecs
    append!(layers,layer(i, x="Signal_E", y="Noise_E", Geom.point, Geom.line))
  end

  p0 = plot(
       layers[1], layers[2],
       Guide.xticks(ticks=ticks_x), Guide.yticks(ticks=ticks_y),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
       Scale.x_log10(labels=x -> @sprintf("%0.0f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.1f", 10^x)),
       Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
       Guide.title(title));
end

function mv_plot(mv_vecs, title="MV", fit=false, save=true)

end

function noise_adjust(a1_df, yc, xc)
    a1_df[:Noise_S_E] = map(x -> sqrt(abs(x^2 - yc^2 - xc^2)), a1_df[:Noise_E]) ;
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~  PTC Examples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
id = ImageDimensions( 100, 100, 100)

# example PTC with normal single output CCD
pix = PixelProperties(0.55, 10, 10, 0)
ac0 = AnalogChain( 0, 1, 10, 7, 1, 0, 19533)
fp0 = FPN(0,0,0,0)
adc0 = ADC(14, 0, 1000000, 2^12-100:2^12+100, 10)

gain = 1/(ac0.cv_gain*ac0.pga_gain*(2^adc0.bit_depth-1)/(adc0.max_volts - adc0.min_volts))
offset = ac0.offset*(2^adc0.bit_depth-1)/(adc0.max_volts - adc0.min_volts)
fwc = Int(round((2^adc0.bit_depth-1) * gain))
q_noise = gain/sqrt(12)

steps = 100
step_size = log10(fwc + 1)/steps

ptc_in, ptc_14 = ptc(pix, ac0, fp0, adc0, id)
ptc_plot([ptc_in, ptc_14])

noise_adjust(ptc_14, 7, q_noise)

p0 = plot(
       layer(ptc_in, x="Signal_E", y="Noise_E", Geom.line, order=1, Theme(default_color=colorant"Red")),
       layer(ptc_14, x="Signal_E", y="Noise_E", Geom.point, Geom.line, order=2, Theme(default_color=colorant"Green")),
       layer(ptc_14, x="Signal_E", y="Noise_S_E", Geom.line, order=3, Theme(default_color=colorant"Blue")),
       Guide.xticks(ticks=[0:1:6;]), Guide.yticks(ticks=[0:0.2:3;]),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
       Scale.x_log10(labels=x -> @sprintf("%0.0f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.1f", 10^x)),
       Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
       Guide.title("Photon Transfer Curve 14-bit, 100K FW"))

 # And saving
 #img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\ptc_14.svg", 6inch, 4inch)
 img = SVG("/Users/rob/Dropbox/Simulator/Graphs/ptc_14_rn_qn.svg", 6inch, 4inch)
 draw(img, p0)
 p0

p0_mv = plot(
        ptc_in, x="Signal_DN", y="Variance_DN", Geom.point, Geom.line, Theme(default_color=colorant"Green"),
        Guide.xticks(ticks=[0:1000:17000;]), Guide.yticks(ticks=[0:200:3000;]),
        style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
        Scale.x_continuous(labels=x -> @sprintf("%0.0f", x)), Scale.y_continuous(labels=x -> @sprintf("%0.0f", x)),
        Coord.cartesian(aspect_ratio=2), #Coord.cartesian(fixed=true),
        Guide.xlabel("Signal(DN)"), Guide.ylabel("Variance(DN^2)"),
        Guide.title("Mean Variance Curve 14-bit, 100K FW"))

p0_mv_zoom = plot(
       ptc_in, x="Signal_DN", y="Variance_DN", Geom.point, Geom.line, Theme(default_color=colorant"Green"),
       Coord.cartesian(xmin=0, xmax=10, ymin=0, ymax=5, aspect_ratio=2),
       Guide.xticks(ticks=[0:2:10;]), Guide.yticks(ticks=[0:1:5;]),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
       Scale.x_continuous(labels=x -> @sprintf("%0.0f", x)), Scale.y_continuous(labels=x -> @sprintf("%0.0f", x)),
       Guide.xlabel("Signal(DN)"), Guide.ylabel("Variance(DN^2)"),
       Guide.title("Mean Variance Curve 14-bit, 100K FW"))

# Regression
mv_fit = fit(LinearModel, @formula(Variance_DN ~ Signal_DN), ptc_14)

# And saving
#img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\mv_zoom_14.svg", 6inch, 4inch)
img = SVG("/Users/rob/Dropbox/Simulator/Graphs/mv_zoom_14.svg", 6inch, 4inch)
draw(img, p0_mv_zoom)
p0_mv_zoom

# example PTC with normal single output EMCCD
a0_em = AnalogChain( 0.04, 128, 10, 80, 1, 0, 19542, 16, 0, 1000000)
steps = 100
gain_em = 1/(a0_em.cv_gain*a0_em.pga_gain*(2^a0_em.bit_depth-1)/(a0_em.max_adc_volts - a0_em.min_adc_volts)) #< system gain w/o EM
fwc_em = Int(round((2^a0_em.bit_depth-1) * gain_em))
step_size_em = log10(fwc_em + 1)/steps
exp, ptc_em_in, ptc_em_16 = ptc(0, 1000, 100, a0_em, id12, pg12)

# slice notation can be used to trim off data points that are not defined or out of range on a log-log axis
p0_em = plot(
       layer(ptc_em_in[1:end-30,:], x="Signal_E", y="Noise_E", Geom.line, order=1, Theme(default_color=colorant"Red")),
       layer(ptc_em_16[2:end-41,:], x="Signal_E", y="Noise_E", Geom.point, Geom.line, order=2, Theme(default_color=colorant"Green")),
       Guide.xticks(ticks=[0:1:4;]), Guide.yticks(ticks=[0:1:2;]),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
       Scale.x_log10(labels=x -> @sprintf("%0.0f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.1f", 10^x)),
       Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
       Guide.title("Photon Transfer Curve EMCCD 16-bit, 100K FW"))

 # And saving
 #img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\ptc_em_16.svg", 6inch, 4inch)
 img = SVG("/Users/rob/Dropbox/Simulator/Graphs/ptc_em_16.svg", 6inch, 4inch)
 draw(img, p0_em)
 p0_em

p0_mv_em = plot(
        ptc_em_in, x="Signal_DN", y="Variance_DN", Geom.point, Geom.line, Theme(default_color=colorant"Green"),
        Guide.xticks(ticks=[0:50:500;]), Guide.yticks(ticks=[0:100:500;]),
        style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
        Scale.x_continuous(labels=x -> @sprintf("%0.0f", x)), Scale.y_continuous(labels=x -> @sprintf("%0.0f", x)),
        Coord.cartesian(aspect_ratio=2), #Coord.cartesian(fixed=true),
        Guide.xlabel("Signal(DN)"), Guide.ylabel("Variance(DN^2)"),
        Guide.title("Mean Variance Curve EMCCD 16-bit, 100K FW"))

p0_mv_zoom_em = plot(
       ptc_in, x="Signal_DN", y="Variance_DN", Geom.point, Geom.line, Theme(default_color=colorant"Green"),
       Coord.cartesian(xmin=0, xmax=10, ymin=0, ymax=5, aspect_ratio=2),
       Guide.xticks(ticks=[0:2:10;]), Guide.yticks(ticks=[0:1:5;]),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=12pt,minor_label_font_size=10pt),
       Scale.x_continuous(labels=x -> @sprintf("%0.0f", x)), Scale.y_continuous(labels=x -> @sprintf("%0.0f", x)),
       Guide.xlabel("Signal(DN)"), Guide.ylabel("Variance(DN^2)"),
       Guide.title("Mean Variance Curve EMCCD 16-bit, 100K FW"))

# Regression
mv_fit = fit(LinearModel, @formula(Variance_DN ~ Signal_DN), ptc_em_in)

# And saving
#img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\ptc1.svg", 6inch, 4inch)
img = SVG("/Users/rob/Dropbox/Simulator/Graphs/mv_em_16.svg", 6inch, 4inch)
draw(img, p0_mv_em)
p0_mv_zoom_em

# example with gain combiner
a1 = AnalogChain( 1, 1, 40, 1.5, 8, 0, 19542)
a2 = AnalogChain( 1, 1, 10, 15, 2, 0, 19542)
adc1 = ADC(11, 0, 4000000,0,0,0)
adc2 = ADC(11, 0, 4000000,0,0,0)
cmb = ADCCombiner(16, 100, 1000, 1100)

# convert offset in uV to DN, system gains in e/DN, FW in electrons
#   - if fractional part is large, offset subtracted results in negative signals
offset_DN = (adc2.adc_offset*(2^adc2.bit_depth-1)/(adc2.max_adc_volts - a2.min_adc_volts))
offset_DN = Int(round(offset_DN))

sg1 = 1/(a1.cv_gain*a1.pga_gain*(2^adc1.bit_depth-1)/(adc1.max_adc_volts - adc1.min_adc_volts))
sg2 = 1/(a2.cv_gain*a2.pga_gain*(2^adc2.bit_depth-1)/(adc2.max_adc_volts - adc2.min_adc_volts))

fw1_electrons = Int(round((2^adc1.bit_depth-1) * sg1))
fw2_electrons = Int(round((2^adc2.bit_depth-1) * sg2))

sgc = fw2_electrons/(2^16-1)

(fw1_electrons, fw2_electrons)

#(min signal(e), max signal(e), steps, analog chain 1, analog chain 2, image dimensions )
ptc_a1, ptc_a2, ptc_c = ptc(0, 21000, 100, a1, a2, id12, pg12)

#pretty plots
rick_morty = Theme(point_size = 3pt, point_size_max = 3pt)
Gadfly.push_theme(rick_morty)

p1 = plot(
       layer(ptc_c, x="Signal_E", y="Noise_E", Geom.point, order=1, Theme(default_color=colorant"Red")),
       layer(ptc_a1, x="Signal_E", y="Noise_E", Geom.line, order=2, Theme(default_color=colorant"Green")),
       layer(ptc_a2, x="Signal_E", y="Noise_E", Geom.line, order=3, Theme(default_color=colorant"Blue")),
       Guide.xticks(ticks=[0:1:6;]), Guide.yticks(ticks=[0:0.2:3;]),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=14pt,minor_label_font_size=12pt),
       Scale.x_log10(labels=x -> @sprintf("%0.0f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.1f", 10^x)),
       Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
       Guide.title("Photon Transfer Curve 14-bit, 100K FW"))

# And saving
#img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\ptc1.svg", 6inch, 4inch)
img = SVG("/Users/rob/Dropbox/Simulator/Graphs/ptc_16_HDR.svg", 6inch, 4inch)
draw(img, p1)
p1

p4 = plot(layer(ptc_c, x="Signal_E", y="Noise_E", Geom.point, order=1, Theme(default_color=colorant"Red")),
           layer(ptc_a1, x="Signal_E", y="Noise_E", Geom.line, order=2, Theme(default_color=colorant"Green")),
           layer(ptc_a2, x="Signal_E", y="Noise_E", Geom.line, order=3, Theme(default_color=colorant"Blue")),
           Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
           Guide.xlabel("Mean Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve: sCMOS HDR"))

p5 = plot(layer(ptc_c, x="Signal_E", y="Variance_E2", Geom.point, order=1, Theme(default_color=colorant"Red")),
           layer(ptc_a1, x="Signal_E", y="Variance_E2", Geom.line, order=2, Theme(default_color=colorant"Green")),
           layer(ptc_a2, x="Signal_E", y="Variance_E2", Geom.line, order=3, Theme(default_color=colorant"Blue")),
           Guide.xlabel("Mean Signal(e)"), Guide.ylabel("Variance(e2)"), Guide.title("Mean Variance Curve: sCMOS HDR"))

p1 = plot(ptc_c, x="Signal_E", y="Noise_E", Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
           Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve CG"))

p2 = plot(ptc_a1, x="Signal_E", y="Noise_E", Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=10),
           Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve HG"))

x_marks = [-1,0,1,2,3]
y_marks = [-1,0,1,2]
p3 = plot(ptc_c, x="Signal_E", y="Noise_E",
       Scale.x_log10(labels=x -> @sprintf("%0.2f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.2f", 10^x)),
       Coord.cartesian(xmin=-1, ymin=-1, ymax=3),
       style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=16pt,minor_label_font_size=14pt),
       Geom.point, Geom.line,
       Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
       Guide.title("Photon Transfer Curve HG"))

p2a = plot(ptc_a2, x="Signal_E", y="Noise_E",
        Scale.x_log10(labels=x -> @sprintf("%0.2f", 10^x)), Scale.y_log10(labels=x -> @sprintf("%0.2f", 10^x)),
        Coord.cartesian(xmin=-1, ymin=-1, ymax=3),
        style(major_label_font="CMU Serif",minor_label_font="CMU Serif", major_label_font_size=16pt,minor_label_font_size=14pt),
        Geom.point, Geom.line,
        Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"),
        Guide.title("Photon Transfer Curve HG"))

p2b = plot(ptc_a2, x="Signal_E", y="Variance_E2", Coord.cartesian, Geom.line, Geom.point,
         Scale.x(minvalue=10), Scale.y(minvalue=10),
         Guide.xlabel("Signal(e)"), Guide.ylabel("Variance(e^2)"), Guide.title("Mean Variance Curve"))

p4 = plot(ptc_c, x="Signal_E", y="Noise_E", Coord.cartesian(xmin=0, ymin=0), Geom.point, Geom.line, order=1,
        Theme(default_color=colorant"Red"), Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=1),
        Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve"))

p5 = plot(ptc_c, x="Signal_E", y="Noise_E", Coord.cartesian(xmin=0, ymin=0), Geom.point, Geom.line, order=1,
        Theme(default_color=colorant"Red"), Scale.x_log10(minvalue=10), Scale.y_log10(minvalue=1),
        Guide.xlabel("Signal(e)"), Guide.ylabel("STD(e)"), Guide.title("Photon Transfer Curve"))

# And saving windows/mac
#img = SVG("C:\\Users\\rlabelle\\Dropbox\\Simulator\\Graphs\\ptc1.svg", 6inch, 4inch)
img = SVG("/Users/rob/Dropbox/Simulator/Graphs/ptc1.svg", 6inch, 4inch)
draw(img, p3)
p3

#display image of a 100e- signal - averaged across 100 frames
em_img = rand(Gamma(10, 100-1+1/10),100,100)
colorview(Gray, em_img/2000)
mean(em_img)
q = DataFrame(ADU = vec(em_img))
plot(q, x="ADU", Geom.histogram)

#~~~~~~~~~~~~~~~~~~~~~~~ importing real image data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
