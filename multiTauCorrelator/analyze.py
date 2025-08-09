import hoomd
from hoomd.custom import Action
from hoomd.logging import log
from hoomd.multiTauCorrelator import _multiTauCorrelator # Import the multi-tau correlator module created by pybind11

# Import the utility function to flatten dictionaries from hoomd.util
# this is needed to handle the logger's nested dictionary structure (see HOOMD documentation of `hoomd.logging.Logger.log()`)
from hoomd.util import _dict_flatten 

class autocorrelate(Action):

    # Indicate that this action will access the pressure tensor
    # (I do not understand fully why, but removing this line
    #  causes the custom action not to work anymore, although
    #  we do not really access the pressure tensor directly in
    #  it)
    flags = [Action.Flags.PRESSURE_TENSOR]

    def __init__(self, quantities: list, logger: hoomd.logging.Logger, filename: str = "autocorrelate.log", eval_period: int = 0, numcorrin: int = 32, p_in: int = 16, m_in: int = 2, normalize: bool = False):
        """
        Initialize the autocorrelator.
        \param quantities: List of quantities to be autocorrelated. They are assumed to be 'HOOMD scalar' quantities.
        \param logger: The logger object containing the quantities to be autocorrelated.
        \param filename: The name of the file to write the correlation data.
        \param period: The frequency at which the autocorrelator is updated.
        \param eval_period: The frequency at which the correlation is evaluated and written to the file.
        \param numcorrin: The number of correlators (levels) in the multi-tau correlator hierarchy.
        \param p_in: The number of time bins for each correlator.
        \param m_in: The decimation factor (or averaging factor) between correlator levels.
                    After every m_in data points at one level, their average is passed to the next level.

        """
        
        super().__init__()
        self.quantities = quantities
        self.logger = logger
        self.filename = filename
        self.eval_period = eval_period
        self.m_corr = None
        self.numcorrin = numcorrin
        self.p_in = p_in
        self.m_in = m_in
        self.firstEvaluation = True
        self.normalize = normalize
    
    @log # placeholder custom logger for reference
    def getFilename(self):
        return self.filename

    def act(self, timestep):
        """
        This is a HOOMD specific method that is called every time the action is triggered.
        This is where the desired action should be coded.
        Other methods can be defined to handle specific actions or calculations and then used inside this method.

        \param timestep: The current simulation timestep.

        This method retrieves the values of the specified quantities from the logger and adds them to the correlator.
        To achieve this, it initializes a separate correlator for each quantity if not already done.

        If eval_period is set, it evaluates the correlator and writes the results to the specified file with the desired frequency.
        If eval_period is 0, it only adds the values to the correlator without writing to the file.

        To write the correlation data to a file, the method `write_to_file` needs to be called
        from the declared `autocorrelate` object in the simulation script. See the function `write_to_file` 
        for the format of the output file.
        """

        if self.m_corr is None:

            self.m_corr = [] # list to hold the correlators for each quantity

            # a dictionary of all the quantities in the logger
            # this is used to check if the quantities are present in the logger
            logger_keys = {key[-1] for key, value in _dict_flatten(self.logger.log()).items()}

            # loop over the quantities in the list and create a correlator for each
            for quantity in self.quantities:
                if quantity not in logger_keys: # check if the quantity is in the logger
                    raise ValueError(f"Quantity '{quantity}' not found in logger keys.")

                # Create a correlator for the quantity (Correlator_Likh object from the C++ module)
                self.m_corr.append(_multiTauCorrelator.Correlator_Likh(self.numcorrin, self.p_in, self.m_in))
                # Initialize the correlator with zero values. This sets up the correlator to be ready for adding values
                self.m_corr[-1].initialize()


        for i, quantity in enumerate(self.quantities):
            # Get the value of the quantity from the logger
            value = self._get_quantity_value(quantity)
            # Add the value to the correlator
            # add() is defined in the C++ module `correlator_likh.cc` and was exposed to Python via pybind11
            self.m_corr[i].add(value)

        if (self.eval_period == 0):
            return

        if (timestep % self.eval_period) == 0:
            self.write_to_file(timestep)


        ### left as a reference for possible future use

        # value = self._get_quantity_value(self.quantities[0])

        # # # print(value, flush=True)
        # # # exit()

        # self.m_corr.add(value)

        # pressure_tensor = self._state._simulation.operations.computes[0].pressure_tensor

        # print(f"Pressure tensor at timestep {timestep}: {pressure_tensor}", flush=True)

        

    def _get_quantity_value(self, quantity):
        """
        Retrieve the value of a specific quantity from the logger.
        \param quantity: The name of the quantity to retrieve.

        Loggers in HOOMD store data in a nested dictionary structure.
        The last dictionary in the nested structure contains the name of the quantity as the key and its value.

        This method flattens the dictionary and retrieves the value for the specified quantity.
        When flattened, the name of the quantity is the last key in the flattened dictionary.
        It returns the first value found for the quantity.
        """

        values = {key[-1]: value[0] for key, value in _dict_flatten(self.logger.log()).items()}
        return values[quantity]


    def write_to_file(self, timestep):
        """
        Write the correlation data to a file.
        \param timestep: The current simulation timestep.

        This method writes loops over the correlators of all the quantities
        and writes the correlation data to the specified file.

        The file format is:
        ```
        correlator evaluated at timestep {timestep}
        timestep,corr_quantity1,corr_quantity2,...
        t0,corr_value1,corr_value2,...
        t1,corr_value1,corr_value2,...
        ...
        ```
        """

        for corr in self.m_corr:
            # Evaluate the correlator to compute the correlation function
            # evaluate() is defined in the C++ module `correlator_likh.cc` and was exposed to Python via pybind11
            corr.evaluate(norm=self.normalize)

        # Open the file in write mode if it's the first evaluation, otherwise append
        if self.firstEvaluation:
            writeMode = "w"
            self.firstEvaluation = False
        else:
            writeMode = "a"

        with open(self.filename, writeMode) as f:

            f.write(f"correlator evaluated at timestep {timestep}\n")

            f.write("timestep")
            for quantity in self.quantities:
                f.write(f",corr_{quantity}")
            f.write("\n")

            # print the number of valid correlation points
            print("number of correlation points:", self.m_corr[0].npcorr, flush=True)

            # Loop over the number of correlation points and
            # write the time and correlation values for each quantity
            for i in range(self.m_corr[0].npcorr):
                f.write(f"{self.m_corr[0].t[i]}")
                for corr in self.m_corr:
                    f.write(f",{corr.f[i]}")
                f.write("\n")
