import warnings

from obspy import UTCDateTime


def write_sacpz(paz, net, sta, cha, fname, include_comments=False):
    """
    From obspy.io.sac.sacpz

    Does not modify normalization_factor, gives same output as converting
    stationxml to dataless using 'stationxml-converter -s', then dataless
    to SAC PZ using 'rdseed -p'

    :param paz: PoleZeroStage representing all PoleZeroStages in the response
    :sta: obspy Station object
    :cha: obspy Channel object
    :fname: filename to write to
    :include_comments: include SAC comments (ruins SEISAN read!)
    """
    sens = cha.response.instrument_sensitivity
    input_unit = sens.input_units.upper()
    if input_unit in ["M", "PA", "PASCALS"]:
        pass
    elif input_unit in ["M/S", "M/SEC"]:
        paz.zeros.append(0j)
    elif input_unit in ["M/S**2", "M/SEC**2"]:
        paz.zeros.extend([0j, 0j])
    else:
        msg = "{}.{}.{}.{} {} ".format(
            net.code,
            sta.code,
            cha.location_code,
            cha.code,
            cha.start_date)
        msg += "has unrecognized input units in "
        msg += "response: {}. Skipping".format(input_unit)
        warnings.warn(msg)
        return
    out = []
    if include_comments is True:
        out = get_comments(paz, net, sta, cha)
    out.append(_paz_to_sacpz_string(paz, sens))
    out.extend(["", ""])
    out = "\n".join(out) + "\n\n"
    with open(fname, "wt") as fh:
        fh.write(out)


def get_comments(paz, net, sta, cha):
    sens = cha.response.instrument_sensitivity
    input_unit = sens.input_units.upper()
    out = []
    now = UTCDateTime()
    out.append("* " + "*" * 50)
    out.append(f"* NETWORK     : {net.code}")
    out.append("* STATION     : %s" % sta.code)
    out.append("* LOCATION    : %s" % cha.location_code)
    out.append("* CHANNEL     : %s" % cha.code)
    out.append("* CREATED     : %s" % now.isoformat())
    out.append("* START       : %s" % cha.start_date.isoformat())
    out.append("* END         : %s" % cha.end_date.isoformat())
    out.append("* DESCRIPTION : %s" % sta.site.name)
    out.append("* LATITUDE    : %s" % (cha.latitude or sta.latitude))
    out.append("* LONGITUDE   : %s" % (cha.longitude or sta.longitude))
    out.append("* ELEVATION   : %s" % (cha.elevation or sta.elevation))
    out.append("* DEPTH       : %s" % cha.depth)
    # DIP in SACPZ served by IRIS SACPZ web service is
    # systematically different from the StationXML entries.
    # It is defined as an incidence angle (like done in SAC for
    # sensor orientation), rather than an actual dip.
    # Add '(SEED)' to clarify that we adhere to SEED convention
    out.append("* DIP (SEED)  : %s" % cha.dip)
    out.append("* AZIMUTH     : %s" % cha.azimuth)
    out.append("* SAMPLE RATE : %s" % cha.sample_rate)
    if input_unit in ['PA', 'PASCALS']:
        out.append(f"* INPUT UNIT  : {input_unit}")
    else:
        out.append("* INPUT UNIT  : M")
    out.append("* OUTPUT UNIT : %s" % sens.output_units)
    out.append("* INSTTYPE    : %s" % paz.description)
    #  out.append("* INSTTYPE    : %s" % cha.sensor.type)
    out.append(f"* INSTGAIN    : {paz.stage_gain:.6e} ({sens.input_units})")
    out.append(f"* SENSITIVITY : {sens.value:.6e} ({sens.input_units})")
    out.append("* A0          : {:.6e} ({})".format(paz.normalization_factor,
                                                    sens.input_units))
    out.append("* " + "*" * 50)
    return out


def _paz_to_sacpz_string(paz, sens):
    """
    Returns SACPZ ASCII text representation of Response.

    COPIED DIRECTLY FROM OBSPY.CORE.INVENTORY.RESPONSE

    :type paz: :class:`PolesZerosResponseStage`
    :param paz: Poles and Zeros response information
    :type instrument_sensitivity: :class:`InstrumentSensitivity`
    :param paz: Overall instrument sensitivity of response
    :rtype: str
    :returns: Textual SACPZ representation of poles and zeros response stage.
    """
    # assemble output string
    out = []
    out.append("ZEROS %i" % len(paz.zeros))
    for c in paz.zeros:
        out.append(" %+.6e %+.6e" % (c.real, c.imag))
    out.append("POLES %i" % len(paz.poles))
    for c in paz.poles:
        out.append(" %+.6e %+.6e" % (c.real, c.imag))
    out.append("CONSTANT {:.6e}".format(paz.normalization_factor * sens.value))
    return "\n".join(out)
