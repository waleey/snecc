class Photon:
    """
    A class representing a photon with its properties and methods.
    """

    def __init__(self, wavelength, pmt, runID, time=0):
        """
        Initialize the Photon with either wavelength or energy.
        :param wavelength: Wavelength of the photon in meters (optional).
        :param energy: Energy of the photon in joules (optional).
        """
        self.wavelength = wavelength  # in meters
        self.pmt = pmt
        self.energy = 1241 / self.wavelength
        self.runID = runID
        self.time = time

    def __repr__(self):
        return f"Photon(wavelength={self.wavelength}, energy={self.energy})"
    
    def update_parameters(self, wavelength=None, pmt = None):
        self.wavelength = wavelength  # in meters
        self.pmt = pmt
        self.energy = 1241 / self.wavelength
    
    def get_energy(self):
        return self.energy
    
    def get_wavelength(self):
        return self.wavelength
    
    def get_pmt(self):
        return self.pmt
    
    def get_runID(self):
        return self.runID
    
    def get_time(self):
        return self.time