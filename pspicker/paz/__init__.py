"""
Pole-Zero representation of an instrument

Follows SEED Pole-Zero protocol for Analog Stages, radian format
(Seed Manual v2.4, Appendix C)

G(f) = S_d * A_0 * prod(s - z_n) / prod(s - p_n) = S_d * A_0 * Hp(s)
where s=j*2*pi*f

The corresponding attributes are:
    S_d: gain factor
    A_0: norm_factor
    G(f): gain at frequency f
    ref_gain = G(f0) = gain at the ref_freq. If there is a passband (a place
                       where the response is flat, ref_freq should be chosen
                       in this band)

The user can specify gain or ref_gain, not both
If norm_factor is not specified, it is automatically calculated so that
    A_0 * abs(Hp(s)) = 1 at f0
as per the SEED recommendation

The user can also change the input_units after creating the class, in which
case the poles, zeros and gain will  be modified if the starting and ending
units are in self.known_units().
"""
from .paz import PAZ
