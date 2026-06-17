from Photon import Photon
import numpy as np
class Event:
    """
    A class describing each neutrino event
    Contains:
    -neutrino energy
    -postrion/electron energy
    -vertex position
    -runID
    -arrival time
    -multiplicity
    -list of photons produced
    """
    def __init__(self, enu, ee, x, y, z, runID, multiplicity = None, time = None, weight = 0.0, min_multiplicity = 1, max_multiplicity = 24):
        self.enu = enu
        self.ee = ee
        self.x = x
        self.y = y
        self.z = z
        self.runID = runID
        self.time = time
        self.multiplicity = multiplicity
        self.photons = []
        self.weight = weight
        self.multiplicity_probability = {
            m: np.zeros(1) 
            for m in range(min_multiplicity, max_multiplicity + 1)
        }

    def __repr__(self):
        return f"Event(enu={self.enu}, elec_energy={self.ee}, position=({self.x}, {self.y}, {self.z}), runID={self.runID}, time={self.time}, multiplicity={self.multiplicity}, n_photons={len(self.photons)})"
    
    def update_parameters(self, enu=None, ee=None, x=None, y=None, z=None, runID=None, multiplicity=None, time=None):
        self.enu = enu if enu is not None else self.enu
        self.ee = ee if ee is not None else self.ee
        self.x = x if x is not None else self.x
        self.y = y if y is not None else self.y
        self.z = z if z is not None else self.z
        self.runID = runID if runID is not None else self.runID
        self.time = time if time is not None else self.time
        self.multiplicity = multiplicity if multiplicity is not None else self.multiplicity

    def add_photon_list(self, photon_list):
        self.photons.extend(photon_list)
    
    def add_photon(self, photon: Photon):
        self.photons.append(photon)

    def calculate_multiplicity(self):
        pmts = []
        for photon in self.photons:
            if photon.get_pmt() not in pmts:
                pmts.append(photon.get_pmt())

        self.multiplicity = len(pmts)

    def setup_multiplicity_probability(self, probabilities):
        for i, p in enumerate(probabilities):
            self.multiplicity_probability[i + 1] = p

    def get_photons(self):
        return self.photons
    
    def get_event_info(self):
        return {
            "enu": self.enu,
            "ee": self.ee,
            "position": (self.x, self.y, self.z),
            "runID": self.runID,
            "time": self.time,
            "multiplicity": self.multiplicity,
            "n_photons": len(self.photons),
            "weight": self.weight
        }
    
    def get_multiplicity(self):
        return self.multiplicity
    
    def get_enu(self):
        return self.enu
    
    def get_ee(self):
        return self.ee
    
    def get_multiplicity_probabilities(self):
        return self.multiplicity_probability
    
    def __setattr__(self, name, value):
        if name == "photons" and not isinstance(value, list):
            raise ValueError("Photons must be a list of Photon objects.")
        super().__setattr__(name, value)