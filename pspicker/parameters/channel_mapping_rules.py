import re


class ChannelMappingRules:
    """
    Rules for mapping channels
    """
    def __init__(self, component_orientation_codes={},
                 write_components_phases={},
                 band_order='GFDCEHSBMLV'):
        """
        :param component_orientation_codes: ComponentOrientationCodes dict
        :param write_components_phases: WriteComponentsPhases dict
        :param band_order: order of preference for band codes (if component
            choice does not yield unique selection)
        """
        self.component_orientation_codes = ComponentOrientationCodes(**component_orientation_codes)
        self.write_components_phases = WriteComponentsPhases(**write_components_phases)
        self.band_order = band_order

    def __str__(self):
        str = "ChannelMappingRules\n"
        str += f"    component_orientation_codes = {self.component_orientation_codes}\n"
        str += f"    write_components_phases = {self.write_components_phases}\n"
        str += f"    band_order = {self.band_order}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return(cls(**thedict))
        
    def seed_id(self, cmp, stream, station):
        """
        Return the seed ID of the trace matching the given station and component

        :param cmp: component to look for (one of 'Z', 'N', 'E', 'H')
        :param stream: Stream containing traces
        :param station: station name
        :returns: None if none found, error if more than one found
        """
        if cmp=='Z':
            cmp_chars = self.component_orientation_codes.Z
        elif cmp=='N':
            cmp_chars = self.component_orientation_codes.N
        elif cmp=='E':
            cmp_chars = self.component_orientation_codes.E
        elif cmp=='H':
            cmp_chars = self.component_orientation_codes.H
        else:
            raise NameError(f"cmp='{cmp}' is not in 'ZNEH'")
        id_match = [tr.get_id() for tr in stream.select(station=station)
                    if re.search(f'[{cmp_chars}]', tr.stats.channel[-1])
                    is not None]
        if len(id_match) > 1:
            for band_code in self.band_order:
                if band_code in [id.split('.')[-1][0] for id in id_match]:
                    id_match = [id for id in id_match
                                if id.split['.'][0] == band_code]
        if len(id_match) == 0:
            return None
        assert len(id_match) == 1,\
            'more than one match found for station = {}, cmp = {}'.format(
                station, cmp_chars)
        return id_match[0]


class ComponentOrientationCodes:
    """
    Component Orientation Codes class
    """
    def __init__(self, Z='Z3', E='E2X', N='N1Y', H='HF'):
        """
        Component characters that will be assigned to a given component
        Each one is provided as a string of possible characters
        
        :param Z: assigned to Z component
        :param E: assigned to E component
        :param N: assigned to N component
        :param H: assigned to H component
        """
        self.Z = Z
        self.E = E
        self.N = N
        self.H = H

    def __str__(self):
        str = "ComponentOrientationCodes\n"
        str += f"    Z = {self.Z}\n"
        str += f"    E = {self.E}\n"
        str += f"    N = {self.N}\n"
        str += f"    H = {self.H}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return(cls(**thedict))

class WriteComponentsPhases:
    """
    Where to write P and S picks [component, phase]
    """
    def __init__(self, S=['N', 'Sg'], P=['Z', 'Pg']):
        """
        First item is the component to write to. Should be
        one of 'ZNEH': output will be written the the component matching
        this value using the ComponentOrientationCodes mappings
        
        Second item is the phase name to give to this arrival 
        
        :param S: S pick parameters  
        :param P: P pick parameters  
        """
        assert S[0] in 'ZNEH', f"S[0]={S[0]}, not in 'ZNEH'"
        assert P[0] in 'ZNEH', f"P[0]={P[0]}, not in 'ZNEH'"
        assert isinstance(S, list)
        assert isinstance(P, list)
        assert len(S) == 2
        assert len(P) == 2
        self.S = S
        self.P = P

    def __str__(self):
        str = "WriteComponentPhases\n"
        str += f"    S = {self.S}\n"
        str += f"    P = {self.P}\n"
        return str

    @classmethod
    def from_dict(cls, thedict):
        return(cls(**thedict))
