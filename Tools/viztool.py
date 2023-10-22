#!/usr/bin/env python3

import matplotlib.ticker as tkr

class FormatScalarFormatter(tkr.ScalarFormatter):
    def __init__(self, fformat='%1.1f', offset=True, mathText=True):
        self.fformat = fformat
        tkr.ScalarFormatter.__init__(self, useOffset=offset, useMathText=mathText)
    def _set_format(self):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format