from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
import PyPR.JSON_Serialization

import numpy as np
import numba
import json

from typing import Iterator, Any, Self
#import system

class FeedbackRegister:
    _seed: np.ndarray[tuple[int],np.dtype[np.uint8]]
    _state: np.ndarray[tuple[int],np.dtype[np.uint8]]
    _prev_state: np.ndarray[tuple[int],np.dtype[np.uint8]]
    fn: FeedbackFunction
    size: int

    #INITIALIATION/DATA:
    def __init__(self,
        seed: int | list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]], 
        fn: FeedbackFunction, 
        compile: bool = False
    ):
        # attributes
        self.fn = fn
        self.size = len(fn)

        # set seed
        self.seed(seed)

        # set states:
        self._state = self._seed.copy()
        self._prev_state = np.zeros_like(self._seed)

        # eat the compilation costs up-front 
        if compile:
            self.fn.compile()
            self.period(limit=1)
            self.clock()
            for state in self.run(1):
                pass
            self.reset()

    def __len__(self) -> int: return self.size

    #REGISTER SEED:
    def seed(self,
        seed: int | list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]]
    ) -> None:
        """Sets the seed for the register.

        The seed value (which is stored in the property `register._seed`) defaults to
        the initial value held in the register. The Register can be set back to that
        state by using `register.reset()`. This method resets that value, using the same
        format as initializing a register.

        We accept either a list or np array to write the bits directly, or an integer, which
        represents the seed in binary. This also enables the use of `random.getrandbits(n)` to
        set an (insecure) random seed of the desired size, which is useful for statistical
        testing.
        
        :param seed: The new value to store as the seed
        :type seed: int | list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]]
        :raises ValueError: If an incompatible seed is passed
        """
        # For a given seed
        if type(seed) == int:
            if seed >= 2**self.size: raise ValueError(f"Seed {seed} larger than register capacity")
            self._seed = np.asarray([int(x) for x in format(seed, f'0{self.size}b')[::-1]], dtype='uint8')
        elif type(seed) == list:
            if len(seed) > self.size: raise ValueError(f"Seed {seed} larger than register capacity")
            self._seed = np.asarray(seed, dtype='uint8')
        elif type(seed) == np.ndarray:
            if len(seed) > self.size: raise ValueError(f"Seed {seed} larger than register capacity")
            self._seed = seed
        else:
            raise ValueError(f'Unexpected seed type {type(seed)}')
        
    def reset(self):
        """Resets the register state (`register._state`) to whatever value is held in `register._seed`
        """
        self._state = self._seed.copy()

    #TYPE CONVERSIONS / CASTING:
    def __str__(self) -> str: return "".join(str(x) for x in self._state[::-1])
    def __int__(self) -> int: return int(str(self),2)
    def __list__(self) -> list: return [int(x) for x in self._state]

    #STATE MANIPULATION:
    def set_state(self,
        state: int | list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]]
    ):
        """Copies the given state into the register's _state array

        We accept either a list or np array to write the bits directly, or an integer, which
        represents the state in binary. This also enables the use of `random.getrandbits(n)` to
        set an random state of the desired size, which is useful for statistical testing.
        
        :param state: The new value to store as the seed
        :type state: int | list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]]
        :raises ValueError: If an incompatible state is passed
        """
        # For a given seed
        if type(state) == int:
            if state >= 2**self.size: raise ValueError(f"State {state} larger than register capacity")
            self._state = np.asarray([int(x) for x in format(state, f'0{self.size}b')[::-1]], dtype='uint8')
        elif type(state) == list:
            if len(state) > self.size: raise ValueError(f"State {state} larger than register capacity")
            self._state = np.asarray(state, dtype='uint8')
        elif type(state) == np.ndarray:
            if len(state) > self.size: raise ValueError(f"State {state} larger than register capacity")
            self._state = state.copy()
        else:
            raise ValueError(f'Unexpected state type {type(state)}')

    def __getitem__(self, key: int) -> int: return self._state[key].copy()
    def __setitem__(self, key: int, val: int): self._state[key] = val
    
    #ITERATION THROUGH REGISTER BITS:
    def __iter__(self): return iter(self._state)
    def __reversed__(self): return reversed(self._state)

    def __copy__(self):
        new_register = FeedbackRegister(0,self.fn.copy())
        new_register._seed = self._seed.copy()
        new_register._state = self._state.copy()
        new_register._prev_state = self._prev_state.copy()
        return new_register
    
    def copy(self):
        new_register = FeedbackRegister(0,self.fn.copy())
        new_register._seed = self._seed.copy()
        new_register._state = self._state.copy()
        new_register._prev_state = self._prev_state.copy()
        return new_register
    
    #CLOCKING AND RUNNING THE REGISTER:
    def clock(self, compiled = True):
        """Update the state held in the state by 1 clock cycle.

        Note that there is both a compiled and uncompiled implementation. Because the compiled version
        is dramatically faster, it is the default, but can be disabled by setting
        `compiled = False`. 

        :param compiled: Whether or not to use the compiled variant.
        :type compiled: bool
        :raises ValueError: If you try to call the compiled implementation while the 
            function has not been compiled
        """
        if compiled and not hasattr(self.fn, "_compiled"):
            raise ValueError(
                "Register Feedback Function is not compiled! " + 
                "Either set compiled = False, or compile the function"
            )
        
        if compiled:
            self._clock_compiled()
        else:
            self._clock_uncompiled()

    def _clock_compiled(self):
        """The compiled branch of the `clock` method.

        As this is the much faster variant, it is the default. When called
        directly it has no input validation, and the responsibility to validate
        that the function is compiled is on the programmer. When routed through
        the generic compile method this is not the case, but when calling this
        branch directly for performance reasons we trust the programmer to know
        what they are doing.
        """
        self._state, self._prev_state = self._prev_state, self._state
        self.fn._compiled_inplace(self._prev_state,self._state)

    def _clock_uncompiled(self):
        """The uncompiled branch of the `clock` method.

        Note that this is dramtically slower than the compiled variant, but may
        be appropriate if you are creating many functions which are only clocked a few
        times.
        """
        # swap pointers to move _state to _prev_state
        self._state, self._prev_state = self._prev_state, self._state

        # overwrite with the new state
        for i in range(len(self.fn.fn_list)):
            self._state[i] = self.fn.fn_list[i].eval(self._prev_state)

    def run(self, limit=None, compiled = True):
        """Returns an iterator which clocks the register for multiple cycles \
            and yields the entire `FeedbackRegister` object.

        Note that there is both a compiled and uncompiled implementation. Because the compiled version
        is dramatically faster, it is the default, but can be disabled by setting
        `compiled = False`. 

        If no limit is passed, the returned iterator will be infinite, and is expected that 
        the user will break the iteration on their end (to avoid an infinite loop).

        If a limit is passed, that many cycles will be returned, counting the initial state
        and excluding the final state. This matches the convention of many other iterators (e.g.
        `range`) to include the start and exclude the end of iteration, and allows seamless connection
        of loops. This is very natural once you get used to it, but worth noting for beginners, as it
        may not be immediately obvious. 
        
        For example, consider `register.run(n)`: this will iterate through the initial state of the 
        register and the next `n-1` states which follow it. If you index the states with the initial
        state being time `0`, then this is states `0` to `n-1`. The registers state after this will be
        the `nth` state, which was not returned in the loop; if another call to `run` is made, this
        current state will be the first one yielded, making multiple sequential calls to `run`
        work together as expected - no states are missed or yielded twice.

        :param limit: the limit on the number of states to iterate through
        :type limit: int | None
        :yield: Yields self on each clock cycle
        :rtype: FeedbackRegister
        :raises ValueError: If you try to call the compiled implementation while the 
            function has not been compiled
        """
        if compiled and not hasattr(self.fn, "_compiled"):
            raise ValueError(
                "Register Feedback Function is not compiled! " + 
                "Either set compiled = False, or compile the function"
            )
        
        if compiled:
            generator = self._run_compiled(limit)
        else:
            generator = self._run_uncompiled(limit)

        for state in generator:
            yield state
        
    def _run_compiled(self, limit: int | None = None) -> Iterator["FeedbackRegister"]:
        """The compiled branch of the `run` method.

        As this is the much faster variant, it is the default. When called
        directly it has no input validation, and the responsibility to validate
        that the function is compiled is on the programmer. When routed through
        the generic compile method this is not the case, but when calling this
        branch directly for performance reasons we trust the programmer to know
        what they are doing. Other than this difference, the behavior is the same as `run`:

        If no limit is passed, the returned iterator will be infinite, and is expected that 
        the user will break the iteration on their end (to avoid an infinite loop).

        If a limit is passed, that many cycles will be returned, counting the initial state
        and excluding the final state. This matches the convention of many other iterators (e.g.
        `range`) to include the start and exclude the end of iteration, and allows seamless connection
        of loops. This is very natural once you get used to it, but worth noting for beginners, as it
        may not be immediately obvious. 
        
        For example, consider `register.run(n)`: this will iterate through the initial state of the 
        register and the next `n-1` states which follow it. If you index the states with the initial
        state being time `0`, then this is states `0` to `n-1`. The registers state after this will be
        the `nth` state, which was not returned in the loop; if another call to `run` is made, this
        current state will be the first one yielded, making multiple calls work together as expected.

        :param limit: the limit on the number of states to iterate through
        :type limit: int | None
        :yield: Yields self on each clock cycle
        :rtype: FeedbackRegister
        """
        # This is here so you can call this branch directly from the object
        # The main implementation is below, outside of the class
        for _ in _run_compiled_(
            self._state,
            self._prev_state,
            self.fn._compiled_inplace,
            limit
        ):
            yield self

    def _run_uncompiled(self, limit: int | None = None) -> Iterator["FeedbackRegister"]:
        """The uncompiled branch of the `run` method.

        Note that this is dramtically slower than the compiled variant, but may
        be appropriate if you are creating many functions which are only clocked a few
        times. Other than this difference, the behavior is the same as `run`:

        If no limit is passed, the returned iterator will be infinite, and is expected that 
        the user will break the iteration on their end (to avoid an infinite loop).

        If a limit is passed, that many cycles will be returned, counting the initial state
        and excluding the final state. This matches the convention of many other iterators (e.g.
        `range`) to include the start and exclude the end of iteration, and allows seamless connection
        of loops. This is very natural once you get used to it, but worth noting for beginners, as it
        may not be immediately obvious. 
        
        For example, consider `register.run(n)`: this will iterate through the initial state of the 
        register and the next `n-1` states which follow it. If you index the states with the initial
        state being time `0`, then this is states `0` to `n-1`. The registers state after this will be
        the `nth` state, which was not returned in the loop; if another call to `run` is made, this
        current state will be the first one yielded, making multiple calls work together as expected.

        :param limit: the limit on the number of states to iterate through
        :type limit: int | None
        :yield: Yields self on each clock cycle
        :rtype: FeedbackRegister
        """
        #number of iterations to run
        if type(limit) == int:
            for _ in range(limit):
                yield self
                self.clock()
                
        #no limit
        elif limit == None:
            while True:
                yield self
                self.clock()


    # PERIOD CALCULATION:
    def period(self, 
        compiled: bool = True, 
        safe: bool = True, 
        limit: int | None = None
    ) -> tuple[int,int] | None:
        """Determine the period of the register with the given state
        
        There are four implementations of this function, depending on the combination
        of the parameters for `safe` and `compiled`, which are explained below. In general,
        unsafe is slightly faster (but less robust) than safe, and compiled is many times
        faster than uncompiled. For general purpose, compiled/safe is the default combination
        which offers a good balance of performance and robustness. The default iteration limit,
        which is used unless an explicit limit is passed, is based on performance to keep the
        runtime roughly under a minute regardless of which implementation is used:
         - For (compiled, unsafe), the default limit is 2^26 cycles.
         - For (compiled, safe), the default limit is 2^25 cycles.
         - For (uncompiled, unsafe), the default limit is 2^16 cycles.
         - For (uncompiled, safe), the default limit is 2^14 cycles.

        :param compiled: A maximum iteration limit,
        :type compiled: bool
        :param safe: if `safe = False`, the period is computed using a naive loop until the state repeats;
            which might fail if the feedback function is non-bijective. if instead `safe = True` (as is default),
            Brent's algorithm for cycle finding is used. Although this adds some overhead, it will never fail.
            This parameter defaults to true, but can be turned off for a performance increase,
            if you are confident the function is bijective
        :type safe: bool
        :param limit: A maximum iteration limit, after which the function will fail and return `None`.
            If `compiled = True`, the limit must be less than 2^64, but since this is impractical to
            check anyway, that is fine. The function exits early upon finfing the period, so if you
            want to have no limit, use a limit such as `2^fn.size`, (or 2^64-1 if the function is 
            larger than 64 bits) which is large enough to catch any cycle. The default values are
            given above.
        :type limit: int, optional
        :return: The period (or None, if the inital state never repeats)
        :rtype: int
        """
        if compiled:
            if not hasattr(self.fn, "_compiled"):
                raise ValueError(
                    "Register Feedback Function is not compiled! " + 
                    "Either set compiled = False, or compile the function"
                )
            if limit and limit >= 2**64-1:
                raise ValueError(
                    "For compiled mode, the largest allowed limit is 2^64-1"
                )
                
        # cases (in order of speed):
        if compiled and not safe:
            if limit == None: limit = 2**26
            return self._period_compiled_unsafe(limit)
        elif compiled and safe: 
            if limit == None: limit = 2**25
            return self._period_compiled_safe(limit)
        elif not compiled and not safe:
            if limit == None: limit = 2**16
            return self._period_uncompiled_unsafe(limit)
        elif not compiled and safe:
            if limit == None: limit = 2**14
            return self._period_uncompiled_safe(limit)
            
    def _period_compiled_unsafe(self, limit: int) -> tuple[int,int] | None:
        """The compiled, unsafe branch of the period function

        This branch is numba compiled, and implements a naive loop, searching for
        repeated states. If the feedback is non-bijective, it may fail and falsly
        report None, even if the limit is much higher than the actual period.

        :param limit: The maximum number of cycles to loop through. The default
            passed in the `period` method is 2^26
        :type limit: int
        :return: Either a tuple containing the period and preperiod (always 0 for this branch),
            or `None` if the method fails.
        :rtype: tuple[int,int] | None
        """

        # This is here so you can call this branch directly from the object
        # The main implementation is below, outside of the class
        return _period_compiled_unsafe_(
            self._state, 
            self.fn._compiled_inplace,
            limit
        )

    def _period_compiled_safe(self,limit: int) -> tuple[int,int] | None:
        """The compiled, safe branch of the period function

        This branch is numba compiled, and implements Brent's cycle-finding algorithm.
        This safely handles non-bijective, feedback, and will simply report the preperiod
        in addition to the period. This is the default branch for `period`.

        :param limit: The maximum number of cycles to loop through. The default
            passed in the `period` method is 2^25
        :type limit: int
        :return: Either a tuple containing the period and preperiod, or `None` if the method fails.
        :rtype: tuple[int,int] | None
        """
        # This is here so you can call  this branch directly from the object
        # The main implementation is below, outside of the class
        return _period_compiled_safe_(
            self._state, 
            self.fn._compiled_inplace,
            limit
        )
  
    def _period_uncompiled_unsafe(self,limit: int) -> tuple[int,int] | None:
        """The uncompiled, unsafe branch of the period function

        This branch is not numba compiled, and so is dramtically slower (but may be
        appropriate if you are creating many small registers in a loop). Like the
        compiled variant, this branch implements a naive loop, searching for
        repeated states. If the feedback is non-bijective, it may fail and falsly
        report None, even if the limit is much higher than the actual period.

        :param limit: The maximum number of cycles to loop through. The default
            passed in the `period` method is 2^16
        :type limit: int
        :return: Either a tuple containing the period and preperiod (always 0 for this branch),
            or `None` if the method fails.
        :rtype: tuple[int,int] | None
        """
        first_state = self._state.copy()
        self.clock()
        count = 1

        while not(all(self._state == first_state)):
            self.clock()
            count += 1

            if count > limit:
                self._state = first_state.copy()
                return None

        self._state = first_state.copy()
        return (count, 0)

    def _period_uncompiled_safe(self,limit: int) -> tuple[int,int] | None:
        """The uncompiled, safe branch of the period function

        This branch is not numba compiled, and so is dramtically slower (but may be
        appropriate if you are creating many small registers in a loop). Like the
        compiled variant, this branch implements Brent's cycle-finding algorithm.
        This safely handles non-bijective, feedback, and will simply report the preperiod
        in addition to the period. This is the default branch for `period`.

        :param limit: The maximum number of cycles to loop through. The default
            passed in the `period` method is 2^25
        :type limit: int
        :return: Either a tuple containing the period and preperiod, or `None` if the method fails.
        :rtype: tuple[int,int] | None
        """
        # Source: slightly modified version of wikipedia's
        # implementation of Brent's algorithm
        slow = FeedbackRegister(self._state.copy(),self.fn)
        fast = FeedbackRegister(self._state.copy(),self.fn)
        fast.clock()

        period = 1
        power = 1

        while not np.all(slow._state == fast._state):
            if power == period:  # time to start a new power of two?
                slow._state = fast._state.copy()
                period = 0
                power *= 2

            fast.clock()
            period += 1

            if period > limit: 
                return None
            
        # Find the position of the first repetition of length lambda
        slow.reset()
        fast.reset()
        for i in range(period):
            fast.clock()
        # The distance between the hare and tortoise is now lambda.
        # Next, the hare and tortoise move at same speed until they agree
        preperiod = 0
        while not np.all(slow._state == fast._state):
            slow.clock()
            fast.clock()
            preperiod += 1
    
        return period, preperiod


    # # JSON / SERIALIZATION
    def generate_ids(self,
        previous_ids: dict[Any, int] | None = None,
        in_place: bool = True
    ) -> dict[Any, int]:
        """Generate a dictionary which maps each object to a unique ID.

        When called with no inputs, this function will create these IDs from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue 
        to add in new IDs, reusing IDs where possible. This means the same function can also be
        used to add onto existing output. For example:
        ```python
        ids = object_1.generate_ids()
        ids = object_2.generate_ids(ids)
        ids = object_3.generate_ids(ids)
        ...
        ```
        For efficiency reasons, this method both mutates the input and returns the mutated
        list of ids by default. This can be disabled by setting the parameter `in_place=False`

        :param previous_ids: the output of previous calls to the function
        :type previous_ids: dict[Any, int] | None
        :param in_place: whether the input list of ids is mutated in place or not. If false, a copy
            is created and returned. If true (by default) the data is modified in place with no copy.
        :type in_place: bool
        :return: A dict which maps each node to a unique id
        :rtype: dict[Any, int]
        """
        if not previous_ids:
            ids = {}
        elif in_place: 
            ids = previous_ids
        else:
            # shallow copy to maintain objects, but new id dict
            ids = {k:v for k,v in previous_ids.items()}

        # create for function:
        ids = self.fn.generate_ids(ids)
        
        # lastly, add id for self:
        ids[self] = max(ids.values()) + 1
        return ids
  
    def _generate_JSON_entry(self,
        ids: dict[Any, int]
    ) -> dict[str, Any]:
        """Create a JSON entry for a given object

        This should return a JSON encoding with all the data necessary to recreate the object
        to an acceptable degree, with a convention matching the corresponding `parse_JSON_entry`.
        This convention is dependent on the class, and this is the pair of methods to overwrite
        to implement JSON serialization for a custom object.

        For `FeedbackRegister` specifically, the `_seed`, `_state`, and `_prev_state` arrays
        are stored as lists, and parse back with the dtype uint8. The `fn` attribute is simply
        stored as a reference to the object's id. The rest of the attributes in `__dict__` are
        copied with no change. 

        :param node_ids: A dictionary which maps each object to a unique id.
        :type node_ids: dict[Any, int]
        :return: A dictionary which represents the JSON data for one object
        :rtype: dict[str, Any]
        """
        # copy class name and non-nested data
        JSON_data = self.__dict__.copy()
        JSON_data['fn'] = ids[self.fn]
        JSON_data['_seed'] = self._seed.tolist()
        JSON_data['_state'] = self._state.tolist()
        JSON_data['_prev_state'] = self._prev_state.tolist()
        return JSON_data

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str,Any],
        parsed_objects: list[Any | None]
    ) -> Self:
        """Parse a JSON entry back into a given object.

        This method is able to parse JSON generated by `_generate_JSON_entry` for the
        corresponding class, according to some convention. This convention may be different,
        and is decided by the class implementer.

        For `FeedbackRegister` specifically, the `_seed`, `_state`, and `_prev_state` arrays
        are stored as lists, and parse back with the dtype uint8. The `fn` attribute is simply
        stored as a reference to the object's id. The rest of the attributes in `__dict__` are
        copied with no change. 

        :param object_data: A dictionary which contains the fields and data of the
            original node object, as generated by `_generate_JSON_entry`.
        :type object_data: dict[str, Any]
        :param parsed_objects: A list which contains the previously parsed objects.
            This can be used to get references to previously stored items, allowing the
            serialization methods to connect the parsed objects together in complex ways.
        :type parsed_objects: list[Any | None]
        :return: The parsed object, with data matching the JSON.
        :rtype: Self
        """
        # intantiate new object:
        new_obj = object.__new__(cls)
        for key,value in object_data.items():
            if key == 'fn':
                new_obj.fn = parsed_objects[value]
            elif key == "_seed":
                new_obj._seed = np.asarray(value,dtype='uint8')
            elif key == "_state":
                new_obj._state = np.asarray(value,dtype='uint8')
            elif key == "_prev_state":
                new_obj._prev_state = np.asarray(value,dtype='uint8')
            else:
                setattr(new_obj,key,value)   
        return new_obj
        
    def to_JSON(self):
        """An alias for `PyPR.JSON_Serialization.generate_JSON(register)`
        
        This can be used in conjunction with `from_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one register. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings on these
        methods for more information on usage and output.

        :return: A JSON object which encodes the input register
        :rtype: dict[str,Any]
        """
        return PyPR.JSON_Serialization.generate_JSON(self)
    
    @classmethod
    def from_JSON(cls, 
        json_object: dict[str,Any]
    ) -> Self:
        """An alias for `PyPR.JSON_Serialization.parse_JSON(json_object)[0]`
        
        This can be used in conjunction with `to_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings on these
        methods for more information on usage and output.

        Although there is no functional difference between `X.from_JSON` and `Y.from_JSON` for two
        classes (`X` and `Y`) which are both serializable, the class you call this method from is used
        to determine type hinting and to clarify the code. Therefore, I choose to throw an error if the
        json encodes a different class than the one you use to decode. This is mostly to enforce 
        readable code and good usage, and to make sure objects are interpreted correctly.

        :param json_object: A dictionary with the expected structure.
        :type json_object: dict[str,Any]
        :return: The FeedbackRegister which was used to create the JSON.
        :rtype: FeedbackRegister
        """
        return_idx = json_object['return order'][0]
        json_class = json_object['objects'][return_idx]['class']
        subclasses = set((
            str(cls)[8:-2] for cls in 
            PyPR.JSON_Serialization.all_subclasses(cls)
        ))
        
        if json_class not in subclasses:
            raise ValueError(
                f"JSON encodes {json_class}, which is not " + 
                f"a subclass of class {str(cls)[8:-2]}"
            )
        
        return PyPR.JSON_Serialization.parse_JSON(json_object)[0]
    
    def to_file(self,
        filename: str
    ) -> None:
        """Writes the output of register.to_JSON to a file with the given filename.
        
        This can be used in conjunction with `from_file` to reduce verbosity and improve readibility
        when you only want to store/parse one object. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings\
        on these methods for more information on usage and output.

        :param filename: A string which will be used as the name of the generated file 
            (must end with the `.json` file extension)
        :type filename: str
        """
        # json files only:
        if filename[-5:] != ".json":
            raise ValueError("Filename must end with the \".json\" file extension")
        
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(cls, 
        filename: str
    ) -> Self:
        """Reads a single register from the file with the given filename.
        
        This can be used in conjunction with `to_file` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        Although there is no functional difference between `X.from_file` and `Y.from_file` for two
        classes (`X` and `Y`) which are both serializable, the class you call this method from is used
        to determine type hinting and to clarify the code. Therefore, I choose to throw an error if the
        json encodes a different class than the one you use to decode. This is mostly to enforce 
        readable code and good usage, and to make sure objects are interpreted correctly.

        :param filename: A string which gives the name of the file to read.
        :type filename: str
        """
        with open(filename, 'r') as f:
            return cls.from_JSON(json.loads(f.read()))


# FULLY COMPILED IMPLEMENTATIONS FOR PERIOD:
u8 = numba.types.uint8
u64 = numba.types.uint64
output_type = numba.types.Optional(numba.types.Tuple([u64,u64]))
update_type = numba.types.FunctionType(numba.void(u8[:],u8[:]))

@numba.njit(output_type(u8[:], update_type, u64))
def _period_compiled_unsafe_(state, update_fn, limit: int) -> tuple[int,int] | None:
    """The compiled, unsafe branch of the period function

    This branch is numba compiled, and implements a naive loop, searching for
    repeated states. If the feedback is non-bijective, it may fail and falsly
    report None, even if the limit is much higher than the actual period.

    :param limit: The maximum number of cycles to loop through. The default
        passed in the `period` method is 2^26
    :type limit: int
    :return: Either a tuple containing the period and preperiod (always 0 for this branch),
        or `None` if the method fails.
    :rtype: tuple[int,int] | None
    """
    prev_state = np.zeros_like(state)
    init_state = state.copy()
    state = state.copy()

    # compiled clock inlined:
    state, prev_state = prev_state, state
    update_fn(prev_state,state)

    for count in range(1,limit+1):
        if np.all(state == init_state):
            return (count,0)
        
        # compiled clock inlined:
        state, prev_state = prev_state, state
        update_fn(prev_state,state)
    return None

@numba.njit(output_type(u8[:], update_type, u64))
def _period_compiled_safe_(state, update_fn, limit: int) -> tuple[int,int] | None:
    """The compiled, safe branch of the period function

    This branch is numba compiled, and implements Brent's cycle-finding algorithm.
    This safely handles non-bijective, feedback, and will simply report the preperiod
    in addition to the period. This is the default branch for `period`.

    :param limit: The maximum number of cycles to loop through. The default
        passed in the `period` method is 2^25
    :type limit: int
    :return: Either a tuple containing the period and preperiod, or `None` if the method fails.
    :rtype: tuple[int,int] | None
    """
    
    # Source: slightly modified version of wikipedia's
    # implementation of Brent's algorithm
    slow_state = state.copy()
    fast_state = state.copy()
    slow_prev_state = np.zeros_like(state)
    fast_prev_state = np.zeros_like(state)

    #fast.clock_compiled() inlined:
    fast_state, fast_prev_state = fast_prev_state, fast_state
    update_fn(fast_prev_state,fast_state)

    period = 1
    power = 1

    while not np.all(slow_state == fast_state):
        if power == period:  # time to start a new power of two?
            slow_state = fast_state.copy()
            period = 0
            power *= 2

        #fast.clock_compiled() inlined:
        fast_state, fast_prev_state = fast_prev_state, fast_state
        update_fn(fast_prev_state,fast_state)
        period += 1

        if period > limit: 
            return None
        
    # Find the position of the first repetition of length λ
    slow_state = state.copy()
    fast_state = state.copy()
    for i in range(period):
        #fast.clock_compiled() inlined:
        fast_state, fast_prev_state = fast_prev_state, fast_state
        update_fn(fast_prev_state,fast_state)
    # The distance between the hare and tortoise is now λ.
    # Next, the hare and tortoise move at same speed until they agree
    preperiod = 0
    while not np.all(slow_state == fast_state):
        #slow.clock_compiled() inlined:
        slow_state, slow_prev_state = slow_prev_state, slow_state
        update_fn(slow_prev_state,slow_state)
        #fast.clock_compiled() inlined:
        fast_state, fast_prev_state = fast_prev_state, fast_state
        update_fn(fast_prev_state,fast_state)
        preperiod += 1

    return period, preperiod

@numba.njit(numba.void(u8[:],u8[:],update_type,numba.types.Optional(u64)))
def _run_compiled_(
    state: np.ndarray,
    prev_state: np.ndarray, 
    update_fn: Any, 
    limit: int | None = None
) -> Iterator[Any]:
    """The compiled branch of the `run` method.

    As this is the much faster variant, it is the default. When called
    directly it has no input validation, and the responsibility to validate
    that the function is compiled is on the programmer. When routed through
    the generic compile method this is not the case, but when calling this
    branch directly for performance reasons we trust the programmer to know
    what they are doing. Other than this difference, the behavior is the same as `run`:

    If no limit is passed, the returned iterator will be infinite, and is expected that 
    the user will break the iteration on their end (to avoid an infinite loop).

    If a limit is passed, that many cycles will be returned, counting the initial state
    and excluding the final state. This matches the convention of many other iterators (e.g.
    `range`) to include the start and exclude the end of iteration, and allows seamless connection
    of loops. This is very natural once you get used to it, but worth noting for beginners, as it
    may not be immediately obvious. 
    
    For example, consider `register.run(n)`: this will iterate through the initial state of the 
    register and the next `n-1` states which follow it. If you index the states with the initial
    state being time `0`, then this is states `0` to `n-1`. The registers state after this will be
    the `nth` state, which was not returned in the loop; if another call to `run` is made, this
    current state will be the first one yielded, making multiple calls work together as expected.

    :param limit: the limit on the number of states to iterate through
    :type limit: int | None
    :yield: Yields self on each clock cycle
    :rtype: FeedbackRegister
    """
    if limit == None:
        # maximum counter allowed - should be more than enough
        for i in range(2**63-1):
            yield

            # XOR swap (can't swap pointers on object so swap values in place)
            for j in range(len(state)):
                state[j] ^= prev_state[j]
                prev_state[j] ^= state[j]
                state[j] ^= prev_state[j]
            update_fn(prev_state,state)

    else:
        for i in range(limit):
            yield

            # XOR swap (can't swap pointers on object so swap values in place)
            for j in range(len(state)):
                state[j] ^= prev_state[j]
                prev_state[j] ^= state[j]
                state[j] ^= prev_state[j]
            update_fn(prev_state,state)
