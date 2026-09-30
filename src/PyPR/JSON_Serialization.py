"""JSON serialization for object graphs that share structure.

A single `generate_JSON(*objs)` call gives every object reachable from its
arguments one integer id, and writes one entry per id. An object referenced
from several places -- a node shared between two functions, the feedback
function a register holds -- is written once and referred to by id everywhere
else, so `parse_JSON` rebuilds the same sharing: every reference comes back as
the same object, not as equal copies.

Classes take part by inheriting `Serializable` and implementing its three
hooks. The mixin supplies the single-object shortcuts (`to_JSON`, `from_JSON`,
`to_file`, `from_file`) and registers each subclass by name, which is how
`parse_JSON` finds the class to rebuild an entry with.
"""
import json
from typing import Any, ClassVar, Self


def qualified_name(cls: type) -> str:
    """The name a class is stored under in the JSON: `module.QualifiedName`.

    :param cls: The class.
    :type cls: type
    :return: The class's module and qualified name, joined by a dot.
    :rtype: str
    """
    return f"{cls.__module__}.{cls.__qualname__}"


class Serializable:
    """Base class for objects stored with `generate_JSON` and `parse_JSON`.

    A subclass implements three hooks, which together fix its JSON encoding:

    - `generate_ids` gives the object, and every serializable object it refers
      to, an id -- reusing the ids already in the map it is handed, which is
      what makes structure shared between objects come back shared.
    - `_generate_JSON_entry` writes the object's own data, with references to
      other serializable objects replaced by their ids.
    - `_parse_JSON_entry` rebuilds the object from that data, resolving ids
      against the objects already rebuilt. Entries are parsed in id order, so
      an object's references must have smaller ids than the object itself;
      `generate_ids` guarantees it by numbering what an object refers to first.

    The default `generate_ids` is correct for an object that refers to no other
    serializable object; a class holding such references overrides it.

    Every subclass is registered under `qualified_name(cls)` when it is
    defined. `parse_JSON` looks classes up in that registry, so the module
    defining a class must be imported before JSON containing it is parsed.
    """

    _registry: ClassVar[dict[str, type]] = {}

    def __init_subclass__(cls, **kwargs: Any):
        super().__init_subclass__(**kwargs)
        Serializable._registry[qualified_name(cls)] = cls

    # -- the contract ------------------------------------------------------
    def generate_ids(self,
        previous_ids: dict[Any, int] | None = None,
        in_place: bool = True
    ) -> dict[Any, int]:
        """Give this object an id, reusing the ids in `previous_ids`.

        Called in sequence, each call continues from the last, so objects shared
        between the calls keep one id::

            ids = object_1.generate_ids()
            ids = object_2.generate_ids(ids)

        This default gives an id to the object alone, which is right for an
        object that refers to no other serializable object. A class holding such
        references overrides it to number them first.

        :param previous_ids: The output of previous calls. Defaults to None.
        :type previous_ids: dict[Any, int] | None
        :param in_place: Whether to add to `previous_ids` itself (the default)
            or to a copy of it.
        :type in_place: bool
        :return: A map from each object to its id.
        :rtype: dict[Any, int]
        """
        if not previous_ids:
            ids = {}
        elif in_place:
            ids = previous_ids
        else:
            ids = dict(previous_ids)

        if self not in ids:
            ids[self] = max(ids.values(), default=-1) + 1
        return ids

    def _generate_JSON_entry(self,
        ids: dict[Any, int]
    ) -> dict[str, Any]:
        """Write the data needed to rebuild this object, in JSON-compatible types.

        References to other serializable objects are written as their ids.

        :param ids: The map from objects to ids made by `generate_ids`.
        :type ids: dict[Any, int]
        :raises NotImplementedError: If not overridden.
        :return: The object's data.
        :rtype: dict[str, Any]
        """
        raise NotImplementedError

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str, Any],
        parsed_objects: list[Any]
    ) -> Self:
        """Rebuild an object from the data `_generate_JSON_entry` wrote.

        :param object_data: The data written for this object.
        :type object_data: dict[str, Any]
        :param parsed_objects: The objects rebuilt so far, indexed by id. Every
            id this object refers to is already filled in.
        :type parsed_objects: list[Any]
        :raises NotImplementedError: If not overridden.
        :return: The rebuilt object.
        :rtype: Self
        """
        raise NotImplementedError

    # -- single-object shortcuts ------------------------------------------
    def to_JSON(self) -> dict[str, Any]:
        """Serialize this object alone: an alias for `generate_JSON(self)`.

        To store several objects whose shared structure should come back shared,
        pass them all to one `generate_JSON` call instead.

        :return: The JSON encoding of this object.
        :rtype: dict[str, Any]
        """
        return generate_JSON(self)

    @classmethod
    def from_JSON(cls,
        json_object: dict[str, Any]
    ) -> Self:
        """Rebuild a single object: an alias for `parse_JSON(json_object)[0]`.

        The class this is called on is checked against the JSON: it must encode
        an instance of this class or a subclass. There is no functional reason
        for the check -- `X.from_JSON` and `Y.from_JSON` would otherwise do the
        same thing -- but it keeps call sites honest about what they load, and
        lets the return type be `Self`.

        :param json_object: JSON produced by `to_JSON` or `generate_JSON`.
        :type json_object: dict[str, Any]
        :raises ValueError: If the JSON encodes a class that is not a subclass
            of this one.
        :return: The first object stored in the JSON.
        :rtype: Self
        """
        return_idx = json_object["return order"][0]
        json_class = _lookup(json_object["objects"][return_idx]["class"])
        if not issubclass(json_class, cls):
            # ValueError, not TypeError: it is what from_JSON has always raised here
            raise ValueError(  # noqa: TRY004
                f"JSON encodes {qualified_name(json_class)}, which is not "
                f"a subclass of {qualified_name(cls)}"
            )
        return parse_JSON(json_object)[0]

    def to_file(self,
        filename: str
    ) -> None:
        """Write `to_JSON()` to a `.json` file.

        :param filename: The file to write; must end in `.json`.
        :type filename: str
        :raises ValueError: If the filename does not end in `.json`.
        """
        if filename[-5:] != ".json":
            raise ValueError("Filename must end with the \".json\" file extension")

        with open(filename, "w") as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(cls,
        filename: str
    ) -> Self:
        """Read a single object from a file written by `to_file`.

        :param filename: The file to read.
        :type filename: str
        :raises ValueError: As for `from_JSON`.
        :return: The object stored in the file.
        :rtype: Self
        """
        with open(filename, "r") as f:
            return cls.from_JSON(json.loads(f.read()))


def _lookup(class_name: str) -> type:
    if class_name not in Serializable._registry:
        raise TypeError(
            f"Type '{class_name}' not recognized: it is not a Serializable class, "
            f"or the module defining it has not been imported"
        )
    return Serializable._registry[class_name]


def generate_JSON(
    *objs: Any
) -> dict[str,Any]:
    """Given a set of objects, generate a JSON file which stores their information.

    Every object reachable from the arguments gets one id, shared structure
    included, and one entry: the class name, and the data its
    `_generate_JSON_entry` writes. A list storing the ids to return helps the
    return of the parse match the structure of the input list.

    :param objs: A variable number of objects to store into the JSON file
    :type objs: Serializable
    :return: The JSON object encoding the data of the input
    :rtype: dict[str, Any]
    """
    # Create the ids:
    obj_ids = {}
    for o in objs:
        # every object needs to specify it's IDs for this to work
        obj_ids = o.generate_ids(obj_ids)

    # create the json:
    num_objs = max(obj_ids.values()) + 1
    obj_entries: list[Any] = [None for i in range(num_objs)]
    for obj, id in obj_ids.items():
        obj_entries[id] = {
            'index': id, # just for human readability
            'class': qualified_name(type(obj)),
            'data': obj._generate_JSON_entry(obj_ids)
        }

    # Annotate with return order:
    return {
        "objects": obj_entries,
        "return order": [obj_ids[obj] for obj in objs]
    }

def parse_JSON(json_object: dict[str,Any]) -> tuple[Any]:
    """Parses a JSON object back into the objects which generated it.

    This method is paired with `generate_JSON`. Entries are rebuilt in id order,
    each by the class registered under its stored name, and an entry's
    references resolve to the objects already rebuilt -- so an object shared in
    the input is a single shared object in the output. The returned tuple
    matches the arguments `generate_JSON` was called with.

    :json_object: A dictionary with the structure generated by `generate_JSON`
    :type: dict[str,Any]
    :raises TypeError: If an entry names a class that is not registered.
    :return: A tuple of objects parsed from the JSON.
    :rtype: tuple[Any]
    """
    # parse object class and data
    return_ids = json_object["return order"]
    json_obj_list = json_object["objects"]
    num_nodes = len(json_obj_list)
    parsed_objects: list[Any] = [None for i in range(num_nodes)]
    for id in range(num_nodes):
        obj_json = json_obj_list[id]
        object_class: Any = _lookup(obj_json['class'])

        # stuff data into new object and add it to the parsed list
        parsed_objects[id] = object_class._parse_JSON_entry(
            obj_json['data'], parsed_objects
        )

    # return the given objects:
    return tuple([parsed_objects[node_id] for node_id in return_ids])
