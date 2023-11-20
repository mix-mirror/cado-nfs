"""
This package is work in progress.

The work plan is:
    1. prime field case (nullspace=right), with arbitrary balancing
        ✅ This works for a few automated tests.
    2. we want to use this to debug the double matrix setting
        - WIP currently
    3. extend to the binary case (nullspace=left). We know that there are
       a few monsters down that path, like m4ri data being hard to
       initialize totally freely from within sage
        ✅ It works in basic cases like mod2_plain.
    4. address the remaining cases like interleaving and so on.
"""

from .BwcParameters import BwcParameters  # noqa: F401
from .BwcMatrix import BwcMatrix  # noqa: F401
from .BwcBalancing import BwcBalancing  # noqa: F401
from .BwcVector import BwcVector  # noqa: F401
from .BwcVectorSet import BwcVectorSet  # noqa: F401
from .BwcXVector import BwcXVector  # noqa: F401
from .BwcCheckData import BwcCheckData  # noqa: F401
from .BwcAFiles import BwcAFiles  # noqa: F401
from .BwcFFiles import BwcFFiles  # noqa: F401
from .BwcSVector import BwcSVector  # noqa: F401
from .BwcSVectorSet import BwcSVectorSet  # noqa: F401
