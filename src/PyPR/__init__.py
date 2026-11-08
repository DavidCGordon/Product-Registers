# PyPR.logging.level = logging.INFO (etc.) sets how much PyPR prints; see
# PyPR.Reporting.LoggingSettings. Modules in the package import the standard
# library's logging as usual -- this name exists only as an attribute here,
# and __all__ leaves it out so `from PyPR import *` cannot shadow a caller's
# own `import logging`.
from PyPR.Reporting import settings as logging
from PyPR.FeedbackRegister import FeedbackRegister

__all__ = ["FeedbackRegister"]
