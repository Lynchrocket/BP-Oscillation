#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True, embedsignature=True
# distutils: language = c++
import numpy as np
import scipy.sparse
from typing import Optional, List, Union
import warnings
import ldpc_modified.helpers.scipy_helpers

cdef BpSparse* Py2BpSparse(pcm):

    cdef int m
    cdef int n
    cdef int nonzero_count = 0

    #check the parity check matrix is the right type
    if isinstance(pcm, np.ndarray) or isinstance(pcm, scipy.sparse.spmatrix):
        pass
    else:
        raise TypeError(f"The input matrix is of an invalid type. Please input\
        a np.ndarray or scipy.sparse.spmatrix object, not {type(pcm)}")

    # Convert to binary sparse matrix and validate input
    pcm = ldpc_modified.helpers.scipy_helpers.convert_to_binary_sparse(pcm)

    # get the parity check dimensions
    m, n = pcm.shape[0], pcm.shape[1]


    # get the number of nonzero entries in the parity check matrix
    if isinstance(pcm,np.ndarray):
        nonzero_count  = int(np.sum( np.count_nonzero(pcm,axis=1) ))
    elif isinstance(pcm,scipy.sparse.spmatrix):
        nonzero_count = int(pcm.nnz)

    # Matrix memory allocation
    cdef BpSparse* cpcm = new BpSparse(m,n,nonzero_count) #creates the C++ sparse matrix object

    #fill sparse matrix
    if isinstance(pcm,np.ndarray):
        for i in range(m):
            for j in range(n):
                if pcm[i,j]==1:
                    cpcm.insert_entry(i,j)
    elif isinstance(pcm,scipy.sparse.spmatrix):
        rows, cols = pcm.nonzero()
        for i in range(len(rows)):
            cpcm.insert_entry(rows[i], cols[i])

    return cpcm

def _native_live_instance_count() -> int:
    return region_adaptive_bp_live_instance_count()



cdef class RegionAdaptiveBpDecoderBase:

    """
    Bp Decoder base class
    """

    def __cinit__(self,pcm, **kwargs):

        error_rate=kwargs.get("error_rate",None)
        error_channel=kwargs.get("error_channel", None)
        max_iter=kwargs.get("max_iter",0)
        bp_method=kwargs.get("bp_method",0)
        ms_scaling_factor=kwargs.get("ms_scaling_factor",1.0)
        schedule=kwargs.get("schedule", 0)
        omp_thread_count = kwargs.get("omp_thread_count", 1)
        random_serial_schedule = kwargs.get("random_serial_schedule", False)
        random_schedule_seed = kwargs.get("random_schedule_seed", 0)
        serial_schedule_order = kwargs.get("serial_schedule_order", None)

        adaptive_max_views = kwargs.get("adaptive_max_views", 4)
        view_thread_count = kwargs.get("view_thread_count", 0)

        assert adaptive_max_views >= 0 and view_thread_count >= 0

        cdef int i, j, nonzero_count
        self.pcm = NULL
        self.bpd = NULL
        self.view_pool = NULL
        self.MEMORY_ALLOCATED = True

        # Matrix memory allocation
        if isinstance(pcm, np.ndarray) or isinstance(pcm, scipy.sparse.spmatrix):
            pass
        else:
            raise TypeError(f"The input matrix is of an invalid type. Please input\
            a np.ndarray or scipy.sparse.spmatrix object, not {type(pcm)}")
        self.pcm = Py2BpSparse(pcm)

        self.m, self.n = pcm.shape[0], pcm.shape[1]

        # allocate vectors for decoder input
        self._error_channel.resize(self.n) #C++ vector for the error channel
        self._syndrome.resize(self.m) #C++ vector for the syndrome
        self._serial_schedule_order = NULL_INT_VECTOR

        if error_channel is not None:
            if len(error_channel) != self.n:
                raise ValueError(f"The error channel vector must have length {self.n}, not {len(error_channel)}.")
            for i in range(self.n):
                assert np.isfinite(error_channel[i]) and 0.0 < error_channel[i] < 1.0
                self._error_channel[i] = error_channel[i]
        elif error_rate is not None:
            if not isinstance(error_rate, float):
                raise ValueError("The `error_rate` parameter must be specified as a single float value.")
            assert np.isfinite(error_rate) and 0.0 < error_rate < 1.0
            for i in range(self.n):
                self._error_channel[i] = error_rate
        else:
            raise ValueError("Please specify the error channel. Either: 1) error_rate: float or 2) error_channel:\
            list of floats of length equal to the block length of the code {self.n}.")

        # Initialise the main decoder and the persistent per-view decoder pool.
        self.bpd = new RegionAdaptiveBpDecoderCpp(self.pcm[0],self._error_channel,0,MINIMUM_SUM,PARALLEL,1.0,1,self._serial_schedule_order,0,False,SYNDROME)
        self.view_pool = new RegionAdaptiveViewDecoderPool()

        ## set the decoder parameters
        self.bp_method = bp_method
        self.max_iter = max_iter
        self.ms_scaling_factor = ms_scaling_factor
        self.schedule = schedule
        self.serial_schedule_order = serial_schedule_order
        self.random_schedule_seed = random_schedule_seed
        self.omp_thread_count = omp_thread_count
        self.random_serial_schedule = random_serial_schedule

        self.adaptive_max_views = adaptive_max_views
        self.view_thread_count = view_thread_count
        self._reset_view_telemetry()

    cdef void _reset_view_telemetry(self):
        self._view_candidates = []
        self._view_iterations = []
        self._view_edge_counts = []
        self._view_validity = []
        self._view_convergence = []

    cdef void _record_view_telemetry(
        self,
        RegionAdaptiveBpDecoderCpp *decoder,
        BpSparse *view_pcm,
    ) except *:
        cdef int i
        cdef bool valid = True
        candidate = np.zeros(self.n, dtype=np.uint8)
        for i in range(self.n):
            candidate[i] = decoder.decoding[i]
        for i in range(self.m):
            if decoder.candidate_syndrome[i] != self._syndrome[i]:
                valid = False
                break
        self._view_candidates.append(candidate)
        self._view_iterations.append(int(decoder.iterations))
        self._view_edge_counts.append(int(view_pcm.entry_count()))
        self._view_validity.append(True if valid else False)
        self._view_convergence.append(True if decoder.converge else False)

    def __dealloc__(self):
        if self.MEMORY_ALLOCATED:
            del self.view_pool
            del self.bpd
            del self.pcm

    @property
    def error_rate(self) -> np.ndarray:
        """
        Returns the current error rate vector.

        Returns:
            np.ndarray: A numpy array containing the current error rate vector.
        """
        out = np.zeros(self.n).astype(float)
        for i in range(self.n):
            out[i] = self.bpd.channel_probabilities[i]
        return out

    @error_rate.setter
    def error_rate(self, value: Optional[float]) -> None:
        """
        Sets the error rate for the decoder.

        Args:
            value (Optional[float]): The error rate value to be set. Must be a single float value.
        """
        if value is not None:
            if not isinstance(value, float):
                raise ValueError("The `error_rate` parameter must be specified as a single float value.")
            assert np.isfinite(value) and 0.0 < value < 1.0
            for i in range(self.n):
                self.bpd.channel_probabilities[i] = value
            self.bpd.refresh_initial_log_prob_ratios()

    @property
    def error_channel(self) -> np.ndarray:
        """
        Returns the current error channel vector.

        Returns:
            np.ndarray: A numpy array containing the current error channel vector.
        """
        out = np.zeros(self.n).astype(float)
        for i in range(self.n):
            out[i] = self.bpd.channel_probabilities[i]
        return out

    @error_channel.setter
    def error_channel(self, value: Union[Optional[List[float]],np.ndarray]) -> None:
        """
        Sets the error channel for the decoder.

        Args:
            value (Optional[List[float]]): The error channel vector to be set. Must have length equal to the block
            length of the code `self.n`.
        """
        if value is not None:
            if len(value) != self.n:
                raise ValueError(f"The error channel vector must have length {self.n}, not {len(value)}.")
            for i in range(self.n):
                assert np.isfinite(value[i]) and 0.0 < value[i] < 1.0
                self.bpd.channel_probabilities[i] = value[i]
            self.bpd.refresh_initial_log_prob_ratios()

    def update_channel_probs(self, value: Union[List[float],np.ndarray]) -> None:
        self.error_channel = value

    @property
    def channel_probs(self) -> np.ndarray:
        out = np.zeros(self.n).astype(float)
        for i in range(self.n):
            out[i] = self.bpd.channel_probabilities[i]
        return out


    @property
    def input_vector_type(self)-> str:
        """
        Returns the current input vector type.

        Returns:
            str: The current input vector type.
        """
        if self.bpd.bp_input_type == SYNDROME:
            return 'syndrome'
        elif self.bpd.bp_input_type == RECEIVED_VECTOR:
            return 'received_vector'
        elif self.bpd.bp_input_type == AUTO:
            return 'auto'
        else:
            raise ValueError(f"The input vector type is invalid. \
                    Please choose from the following methods: \
                    'input_vector_type=syndrome', 'input_vector_type=received_vector'")


    @input_vector_type.setter
    def input_vector_type(self, input_type: str):
        """
        Sets the input vector type.

        Args:
            input_type (str): The input vector type to be set. Must be either 'syndrome' or 'received_vector'.
        """
        if input_type.lower() in ['auto', 'a', '2']:
            if self.m == self.n:
                raise ValueError("Please specify the input vector type. Either: 1) input_vector_type: 'syndrome' or 2) input_vector_type:\
                'received_vector'.")
            else:
                self.bpd.bp_input_type = AUTO

        elif input_type.lower() in ['syndrome', 's', '0']:
            self.bpd.bp_input_type = SYNDROME
        elif input_type.lower() in ['received_vector', 'r', '1']:
            self.bpd.bp_input_type = RECEIVED_VECTOR
        else:
            raise ValueError(f"The input vector type '{input_type}' is invalid. \
                    Please choose from the following methods: \
                    'input_vector_type=syndrome', 'input_vector_type=received_vector'")


    @property
    def log_prob_ratios(self) -> np.ndarray:
        """
        Returns the current log probability ratio vector.

        Returns:
            np.ndarray: A numpy array containing the current log probability ratio vector.
        """
        out = np.zeros(self.n)
        for i in range(self.n):
            out[i] = self.bpd.log_prob_ratios[i]
        return out

    @property
    def auxiliary_log_prob_ratios(self) -> np.ndarray:
        cdef RegionAdaptiveBpDecoderCpp *decoder
        out = np.zeros(self.n)
        if self.view_pool.active_view_count() == 0:
            for i in range(self.n):
                out[i] = self.bpd.initial_log_prob_ratios[i]
            return out
        decoder = &self.view_pool.decoder(self.view_pool.active_view_count() - 1)
        for i in range(self.n):
            out[i] = decoder.log_prob_ratios[i]
        return out

    @property
    def converge(self) -> bool:
        """
        Returns whether the decoder has converged or not.

        Returns:
            bool: True if the decoder has converged, False otherwise.
        """
        return self.bpd.converge

    @property
    def auxiliary_converge(self) -> bool:
        if self.view_pool.active_view_count() == 0:
            return False
        return self.view_pool.decoder(self.view_pool.active_view_count() - 1).converge

    @property
    def iter(self) -> int:
        """
        Returns the number of iterations performed by the decoder.

        Returns:
            int: The number of iterations performed by the decoder.
        """
        return self.bpd.iterations


    @property
    def check_count(self) -> int:
        """
        Returns the number of rows of the parity check matrix.

        Returns:
            int: The number of rows of the parity check matrix.
        """
        return self.bpd.pcm.m

    @property
    def bit_count(self) -> int:
        """
        Returns the number of columns of the parity check matrix.

        Returns:
            int: The number of columns of the parity check matrix.
        """
        return self.bpd.pcm.n

    @property
    def max_iter(self) -> int:
        """
        Returns the maximum number of iterations allowed by the decoder.

        Returns:
            int: The maximum number of iterations allowed by the decoder.
        """
        return self.bpd.maximum_iterations

    @max_iter.setter
    def max_iter(self, value: int) -> None:
        """
        Sets the maximum number of iterations allowed by the decoder.

        Args:
            value (int): The maximum number of iterations allowed by the decoder.

        Raises:
            ValueError: If value is not a positive integer.
        """
        if not isinstance(value, int):
            raise ValueError("max_iter input parameter is invalid. This must be specified as a positive int.")
        if value < 0:
            raise ValueError(f"max_iter input parameter must be a positive int. Not {value}.")
        self.bpd.maximum_iterations = value if value != 0 else self.n

    @property
    def bp_method(self) -> str:
        """
        Returns the belief propagation method used.

        Returns:
            str: The belief propagation method used. Possible values are 'product_sum' or 'minimum_sum'.
        """
        if self.bpd.bp_method == PRODUCT_SUM:
            return 'product_sum'
        elif self.bpd.bp_method == MINIMUM_SUM:
            return 'minimum_sum'
        else:
            raise ValueError(f"BP method is invalid. \
                    Please choose from the following methods: \
                    'product_sum', 'minimum_sum'")

    @bp_method.setter
    def bp_method(self, value: Union[str,int]) -> None:
        """
        Sets the belief propagation method used.

        Args:
            value (str): The belief propagation method to use. Possible values are 'product_sum' or 'minimum_sum'.

        Raises:
            ValueError: If value is not a valid option.
        """
        if str(value).lower() in ['prod_sum', 'product_sum', 'ps', '0', 'prod sum']:
            self.bpd.bp_method = PRODUCT_SUM
        elif str(value).lower() in ['min_sum', 'minimum_sum', 'ms', '1', 'minimum sum', 'min sum']:
            self.bpd.bp_method = MINIMUM_SUM
        else:
            raise ValueError(f"BP method '{value}' is invalid. \
                    Please choose from the following methods: \
                    'product_sum', 'minimum_sum'")

    @property
    def schedule(self) -> str:
        """
        Returns the scheduling method used.

        Returns:
            str: The scheduling method used. Possible values are 'parallel' or 'serial'.
        """
        if self.bpd.schedule == PARALLEL:
            return 'parallel'
        elif self.bpd.schedule == SERIAL:
            return 'serial'
        elif self.bpd.schedule == SERIAL_RELATIVE:
            return 'serial_relative'
        else:
            raise ValueError(f"The BP schedule method is invalid. \
                    Please choose from the following methods: \
                    'schedule=parallel', 'schedule=serial', 'schedule=serial_relative'")

    @schedule.setter
    def schedule(self, value: Union[str,int]) -> None:
        """
        Sets the scheduling method used.

        Args:
            value (str): The scheduling method to use. Possible values are 'parallel' or 'serial'.

        Raises:
            ValueError: If value is not a valid option.
        """
        if str(value).lower() in ['parallel','p','0']:
            self.bpd.schedule = PARALLEL
        elif str(value).lower() in ['serial','s','1']:
            self.bpd.schedule = SERIAL
        elif str(value).lower() in ['serial_relative', 'sr', '2']:
            self.bpd.schedule = SERIAL_RELATIVE
        else:
            raise ValueError(f"The BP schedule method '{value}' is invalid. \
                    Please choose from the following methods: \
                    'schedule=parallel', 'schedule=serial', 'schedule=serial_relative'")

    @property
    def serial_schedule_order(self) -> Union[None, np.ndarray]:
        """
        Returns the serial schedule order.

        Returns:
            Union[None, np.ndarray]: The serial schedule order as a numpy array, or None if no schedule has been set.
        """
        if self.bpd.serial_schedule_order.size() == 0:
            return None

        out = np.zeros(self.n).astype(int)
        for i in range(self.n):
            out[i] = self.bpd.serial_schedule_order[i]
        return out

    @serial_schedule_order.setter
    def serial_schedule_order(self, value: Union[None, List[int], np.ndarray]) -> None:
        """
        Sets the serial schedule order.

        Args:
            value (Union[None, List[int]]): The serial schedule order to set. Must have length equal to the block
            length of the code `self.n`.

        Raises:
            Exception: If value does not have the correct length.
            ValueError: If value contains an invalid integer value.
        """
        if value is None:
            self._serial_schedule_order = NULL_INT_VECTOR
            for i in range(self.n):
                self.bpd.serial_schedule_order[i] = i
            return
        if not len(value) == self.n:
            raise Exception("Input error. The `serial_schedule_order` input parameter must have length equal to the length of the code.")
        for i in range(self.n):
            if not isinstance(value[i], (int, np.int64, np.int32)) or value[i] < 0 or value[i] >= self.n:
                raise ValueError(f"serial_schedule_order[{i}] is invalid. It must be a non-negative integer less than {self.n}.")
        if len({int(item) for item in value}) != self.n:
            raise ValueError("serial_schedule_order must be a permutation of bit indices.")
        for i in range(self.n):
            self.bpd.serial_schedule_order[i] = value[i]
        self.random_serial_schedule = False

    @property
    def ms_scaling_factor(self) -> float:
        """Get the scaling factor for minimum sum method.

        Returns:
            float: The current scaling factor.
        """
        return self.bpd.ms_scaling_factor

    @ms_scaling_factor.setter
    def ms_scaling_factor(self, value: float) -> None:
        """Set the scaling factor for minimum sum method.

        Args:
            value (float): The new scaling factor.

        Raises:
            TypeError: If the input value is not a float.
        """
        if not isinstance(value, (float,int)):
            raise TypeError("The ms_scaling factor must be specified as a float")
        self.bpd.ms_scaling_factor = value

    @property
    def omp_thread_count(self) -> int:
        """Get the number of OpenMP threads.

        Returns:
            int: The number of threads used.
        """
        if self.bpd.omp_thread_count != 1:
            warnings.warn("The OpenMP functionality is not yet implemented")
        return self.bpd.omp_thread_count

    @omp_thread_count.setter
    def omp_thread_count(self, value: int) -> None:
        """Set the number of OpenMP threads.

        Args:
            value (int): The number of threads to use.

        Raises:
            TypeError: If the input value is not an integer or is less than 1.
        """
        if not isinstance(value, int) or value < 1:
            raise TypeError("The omp_thread_count must be specified as a\
            positive integer.")
        self.bpd.set_omp_thread_count(value)
        if self.bpd.omp_thread_count != 1:
            warnings.warn("The OpenMP functionality is not yet implemented")

    @property
    def random_schedule_seed(self) -> int:
        """Get the value of random_schedule_seed.

        Returns:
            int: The current value of random_schedule_seed.
        """
        return self.bpd.random_schedule_seed

    @random_schedule_seed.setter
    def random_schedule_seed(self, value: int) -> None:
        """Set the value of random_schedule_seed.

        Args:
            value (int): The new value of random_schedule_seed.

        Raises:
            ValueError: If the input value is not a postive integer.
        """
        if not isinstance(value, int) or value < -2:
            raise ValueError("The value of random_schedule_seed must\
            be a positive integer. Set as -1 to disable to the random\
            schedule. Set as 0 to use the system clock.")

        self.bpd.random_serial_schedule = True
        self.bpd.set_random_schedule_seed(value)

    @property
    def random_serial_schedule(self) -> bool:
        """
        Returns whether the random serial schedule is enabled.

        Returns:
            bool: True if random serial schedule is enabled, False otherwise.
        """
        return self.bpd.random_serial_schedule

    @random_serial_schedule.setter
    def random_serial_schedule(self, value: bool) -> None:
        """
        Sets whether the random serial schedule is enabled.

        Args:
            value (int): True to enable random serial schedule, False to disable it.

        Raises:
            ValueError: If random serial schedule is enabled while a fixed serial schedule is set.
        """
        # if not isinstance(value, bool):
        #     raise ValueError("The random_serial_schedule must be a boolean value.")
        self.bpd.random_serial_schedule = value

    @property
    def sign_flip_count(self) -> np.ndarray:
        """
        Returns the sign flip count vector.

        Returns:
            np.ndarray: A numpy array containing the sign flip count for each bit.
        """
        out = np.zeros(self.n).astype(int)
        for i in range(self.n):
            out[i] = self.bpd.sign_flip_count[i]
        return out

    @property
    def llr_range(self) -> np.ndarray:
        """
        Returns the LLR range vector.

        Returns:
            np.ndarray: A numpy array containing the LLR range for each bit.
        """
        out = np.zeros(self.n).astype(float)
        for i in range(self.n):
            out[i] = self.bpd.llr_max[i] - self.bpd.llr_min[i]
        return out

cdef class RegionAdaptiveBpDecoder(RegionAdaptiveBpDecoderBase):
    """
    Belief propagation decoder for binary linear codes.

    This class provides an implementation of belief propagation decoding for binary linear codes. The decoder uses a sparse
    parity check matrix to decode received codewords. The decoding algorithm can be configured using various parameters,
    such as the belief propagation method used, the scheduling method used, and the maximum number of iterations.

    Parameters
    ----------
    pcm : Union[np.ndarray, spmatrix]
        The parity check matrix of the binary linear code, represented as a NumPy array or a SciPy sparse matrix.
    error_rate : Optional[float], optional
        The initial error rate for the decoder, by default None.
    error_channel : Optional[List[float]], optional
        The initial error channel probabilities for the decoder, by default None.
    max_iter : Optional[int], optional
        The maximum number of iterations allowed for decoding, by default 0 (adaptive).
    bp_method : Optional[str], optional
        The belief propagation method to use: 'product_sum' or 'minimum_sum', by default 'minimum_sum'.
    ms_scaling_factor : Optional[float], optional
        The scaling factor for the minimum sum method, by default 1.0.
    schedule : Optional[str], optional
        The scheduling method for belief propagation: 'parallel', 'serial', or 'serial_relative'. By default 'parallel'.
    omp_thread_count : Optional[int], optional
        The number of OpenMP threads to use, by default 1.
    random_schedule_seed : Optional[int], optional
        The seed for the random serial schedule, by default 0. If set to 0, the seed is set according to the system clock.
    serial_schedule_order : Optional[List[int]], optional
        The custom order for serial scheduling, by default None.
    random_serial_schedule : bool, optional
        Whether to enable random serial scheduling. If True, the serial schedule order is randomized in each iteration.
        By default False.
    input_vector_type: str, optional
        Use this parameter to specify the input type. Choose either: 1) 'syndrome' or 2) 'received_vector' or 3) 'auto'.
        Note, it is only necessary to specify this value when the parity check matrix is square. When the
        parity matrix is non-square, the input vector type is inferred automatically from its length.
    """

    def __cinit__(self, pcm: Union[np.ndarray, scipy.sparse.spmatrix], error_rate: Optional[float] = None,
                 error_channel: Optional[Union[np.ndarray,List[float]]] = None, max_iter: Optional[int] = 0, bp_method: Optional[str] = 'minimum_sum',
                 ms_scaling_factor: Optional[Union[float,int]] = 1.0, schedule: Optional[str] = 'parallel', omp_thread_count: Optional[int] = 1,
                 random_schedule_seed: Optional[int] = 0, serial_schedule_order: Optional[List[int]] = None, input_vector_type: str = "auto", random_serial_schedule: bool = False,
                 adaptive_max_views: int = 4, view_thread_count: int = 0, **kwargs):

        # for key in kwargs.keys():
        #     if key not in ["channel_probs"]:
        #         raise ValueError(f"Unknown parameter '{key}' passed to the RegionAdaptiveBpDecoder constructor.")

        self.input_vector_type = input_vector_type
        self._received_vector.resize(self.n) #C++ vector for the received vector

        pass

    def __init__(self, pcm: Union[np.ndarray, scipy.sparse.spmatrix], error_rate: Optional[float] = None,
                 error_channel: Optional[Union[np.ndarray,List[float]]] = None, max_iter: Optional[int] = 0, bp_method: Optional[str] = 'minimum_sum',
                 ms_scaling_factor: Optional[Union[float,int]] = 1.0, schedule: Optional[str] = 'parallel', omp_thread_count: Optional[int] = 1,
                 random_schedule_seed: Optional[int] = 0, serial_schedule_order: Optional[List[int]] = None,
                 input_vector_type: str = "auto", random_serial_schedule: bool = False,
                 adaptive_max_views: int = 4, view_thread_count: int = 0, **kwargs):

        pass

    def decode(self, input_vector: np.ndarray) -> np.ndarray:
        """Decode with one main BP run and complementary auxiliary BP views."""

        self._reset_view_telemetry()

        if self.bpd.bp_input_type == SYNDROME and len(input_vector) != self.m:
            raise ValueError(f"The input_vector must have length {self.m} (for syndrome decoding). Not length {len(input_vector)}.")
        elif self.bpd.bp_input_type == RECEIVED_VECTOR and len(input_vector) != self.n:
            raise ValueError(f"The input_vector must have length {self.n} (for received vector decoding). Not length {len(input_vector)}.")
        elif self.bpd.bp_input_type == AUTO and len(input_vector) not in (self.m, self.n):
            raise ValueError(f"The input_vector must have length {self.m} (for syndrome decoding) or length {self.n} (for received vector decoding). Not length {len(input_vector)}.")

        cdef int i
        cdef int view_index
        cdef bool zero_input_vector = True
        cdef RegionAdaptiveBpDecoderCpp *view_decoder
        cdef BpSparse *view_pcm
        DTYPE = input_vector.dtype
        cdef int len_input_vector = len(input_vector)

        if self.bpd.bp_input_type == SYNDROME or (self.bpd.bp_input_type == AUTO and len_input_vector == self.m):
            for i in range(len_input_vector):
                self._syndrome[i] = input_vector[i]
                if self._syndrome[i]: zero_input_vector = False
            if zero_input_vector:
                self.bpd.reset_zero_syndrome_state()
                self.view_pool.reset()
                self.bpd.reset_adaptive_multi_view_state()
                self.bpd.consider_adaptive_candidate(
                    self._syndrome,
                    self.bpd.decoding,
                    self.bpd.candidate_syndrome,
                    0,
                )
                return np.zeros(self.bit_count, dtype=DTYPE)

            self.bpd.decode(self._syndrome)
            self.bpd.prepare_region_redundant_views(
                self._syndrome,
                self.adaptive_max_views,
            )
            self.bpd.adaptive_main_bp_run_count = 1
            self.view_pool.reset()
            self.bpd.consider_adaptive_candidate(
                self._syndrome,
                self.bpd.decoding,
                self.bpd.candidate_syndrome,
                0,
            )
            self._record_view_telemetry(self.bpd, self.pcm)

            self.view_pool.run(self.bpd[0], self._syndrome, self.view_thread_count)
            for view_index in range(self.view_pool.active_view_count()):
                view_decoder = &self.view_pool.decoder(view_index)
                view_pcm = &self.view_pool.pcm(view_index)
                self._record_view_telemetry(view_decoder, view_pcm)
                self.bpd.consider_adaptive_candidate(
                    self._syndrome,
                    view_decoder.decoding,
                    view_decoder.candidate_syndrome,
                    view_index + 1,
                )

        else:
            for i in range(len_input_vector):
                self._received_vector[i] = input_vector[i]
                if self._received_vector[i]: zero_input_vector = False
            if zero_input_vector:
                self.bpd.reset_zero_syndrome_state()
                self.view_pool.reset()
                self.bpd.reset_adaptive_multi_view_state()
                return np.zeros(self.bit_count, dtype=DTYPE)
            self.bpd.decode(self._received_vector)
            self.view_pool.reset()
            self.bpd.reset_adaptive_multi_view_state()
            self.bpd.adaptive_main_bp_run_count = 1
            self.bpd.adaptive_candidate_count = 1
            self.bpd.adaptive_selected_view = 0
            self.bpd.adaptive_final_decoding = self.bpd.decoding

        out = np.zeros(self.n, dtype=DTYPE)
        for i in range(self.n):
            out[i] = self.bpd.adaptive_final_decoding[i]
        return out


    @property
    def decoding(self) -> np.ndarray:
        """Return the candidate selected across all BP views."""
        out = np.zeros(self.n, dtype=int)
        for i in range(self.n):
            out[i] = self.bpd.adaptive_final_decoding[i]
        return out

    @property
    def main_decoding(self) -> np.ndarray:
        out = np.zeros(self.n, dtype=int)
        for i in range(self.n):
            out[i] = self.bpd.decoding[i]
        return out

    @property
    def auxiliary_decoding(self) -> np.ndarray:
        cdef RegionAdaptiveBpDecoderCpp *decoder
        out = np.zeros(self.n, dtype=int)
        if self.view_pool.active_view_count() == 0:
            return out
        decoder = &self.view_pool.decoder(self.view_pool.active_view_count() - 1)
        for i in range(self.n):
            out[i] = decoder.decoding[i]
        return out

    @property
    def instability_scores(self) -> np.ndarray:
        self.bpd.ensure_instability_scores()
        out = np.zeros(self.n, dtype=float)
        for i in range(self.n):
            out[i] = self.bpd.instability_scores[i]
        return out

    @property
    def view_candidates(self) -> np.ndarray:
        if not self._view_candidates:
            return np.zeros((0, self.n), dtype=np.uint8)
        return np.asarray(self._view_candidates, dtype=np.uint8).copy()

    @property
    def view_iterations(self) -> np.ndarray:
        return np.asarray(self._view_iterations, dtype=np.int64)

    @property
    def view_edge_counts(self) -> np.ndarray:
        return np.asarray(self._view_edge_counts, dtype=np.int64)

    @property
    def view_validity(self) -> np.ndarray:
        return np.asarray(self._view_validity, dtype=np.bool_)

    @property
    def view_convergence(self) -> np.ndarray:
        return np.asarray(self._view_convergence, dtype=np.bool_)

    @property
    def region_stats(self) -> dict:
        return {
            "prepared_view_count": int(self.bpd.adaptive_prepared_view_count()),
            "executed_view_count": int(self.bpd.adaptive_executed_view_count),
            "selected_rows_per_view": [
                int(self.bpd.adaptive_selected_row_count(i))
                for i in range(self.bpd.adaptive_prepared_view_count())
            ],
            "automatic_weight_limit_per_view": [
                int(self.bpd.adaptive_weight_limit(i))
                for i in range(self.bpd.adaptive_prepared_view_count())
            ],
            "selected_row_weights_per_view": [
                [
                    int(self.bpd.adaptive_selected_row_weight(i, j))
                    for j in range(self.bpd.adaptive_selected_row_count(i))
                ]
                for i in range(self.bpd.adaptive_prepared_view_count())
            ],
            "selected_row_actionability_per_view": [
                [
                    float(self.bpd.region_selected_row_actionability(i, j))
                    for j in range(self.bpd.adaptive_selected_row_count(i))
                ]
                for i in range(self.bpd.adaptive_prepared_view_count())
            ],
            "selected_row_cycle_density_per_view": [
                [
                    float(self.bpd.region_selected_row_cycle_density(i, j))
                    for j in range(self.bpd.adaptive_selected_row_count(i))
                ]
                for i in range(self.bpd.adaptive_prepared_view_count())
            ],
            "core_variable_count": int(self.bpd.region_core_variable_count()),
            "conflict_check_count": int(self.bpd.region_conflict_check_count()),
            "anchor_variable_count": int(self.bpd.region_anchor_variable_count()),
            "unstable_variable_count": int(self.bpd.adaptive_unstable_variable_count),
            "feasible_unstable_variable_count": int(self.bpd.adaptive_feasible_unstable_variable_count),
            "covered_unstable_variable_count": int(self.bpd.adaptive_covered_unstable_variable_count),
            "candidate_pair_count": int(self.bpd.adaptive_candidate_pair_count()),
            "candidate_count": int(self.bpd.adaptive_candidate_count),
            "valid_candidate_count": int(self.bpd.adaptive_valid_candidate_count),
            "selected_view": int(self.bpd.adaptive_selected_view),
            "best_candidate_cost": float(self.bpd.adaptive_best_candidate_cost),
            "fallback_residual_weight": int(self.bpd.adaptive_fallback_residual_weight),
            "selected_candidate_valid": True if self.bpd.adaptive_selected_candidate_valid else False,
            "view_worker_count": int(self.view_pool.last_worker_count()),
            "bp_run_count": int(
                self.bpd.adaptive_main_bp_run_count + self.bpd.adaptive_executed_view_count
            ),
        }
