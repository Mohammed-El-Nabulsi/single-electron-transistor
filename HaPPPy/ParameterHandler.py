

class ParameterHandler():
    """A simple class to store parameters needed by more than one project."""

    temperature = None

    def setTemperature(self, temperature):
        """Sets the temperature."""
        self.temperature = temperature

    def getTemperature(self):
        """Returns the stored temperature."""
        return self.temperature
        
