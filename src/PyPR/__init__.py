# PyPR.logging.level = logging.INFO (etc.) sets how much PyPR prints; see
# PyPR.Reporting.LoggingSettings. Modules in the package import the standard
# library's logging as usual -- this name exists only as an attribute here.
from PyPR.Reporting import settings as logging
from PyPR.FeedbackRegister import FeedbackRegister
