#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Ionic enumerations** (i.e., :class:`enum.Enum` subclasses required to
differentiate between various types of ions).
'''

# ....................{ IMPORTS                           }....................
from betse.util.type import enums
from betse.util.type.types import GeneratorType

# ....................{ ITERATORS                         }....................
def iter_ion_names() -> GeneratorType:
    '''
    Generator yielding the unmodified name of each supported type of ion
    specified by the :class:`IonType` enumeration (in lexicographic order).

    Returns
    ----------
    GeneratorType
        Generator yielding each such name.
    '''

    # Let there be ion.
    yield from enums.iter_names(IonType)

# ....................{ ENUMS                             }....................
IonType = enums.make_enum(
    class_name='IonType',
    member_names=('Na', 'K', 'Ca', 'Cl', 'M', 'P',),
    is_ordered=True,
    doc='''
Ordered enumeration of all supported types of **ions** (i.e., atoms or
molecules having lost or gained electrons and hence net electric charge).

For brevity, the names of these types are intentionally condensed to their
standard abbreviations rather than elongated to their human-readable nouns
(e.g., ``Na`` rather than ``SODIUM``). Doing so also improves interoperability
with the remainder of the codebase, which first invented this nomenclature.

Attributes
----------
Na : enum
    Sodium cation Na+.
K : enum
    Potassium cation K+.
Ca : enum
    Calcium cation Ca2+.
Cl : enum
    Chloride anion Cl-.
M : enum
    Unidentified anion M-, synthetically manufactured to enforce
    charge-balancing across the environmental matrix for both simulation
    stability and correctness.
P : enum
    Anionic protein P-.
''')
