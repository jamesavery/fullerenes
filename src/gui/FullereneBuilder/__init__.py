# (c) Copyright QuantumWise A/S
# See the file "LICENSE" in the AddOns directory of the ATK distribution
# for the full license governing this code.
import FullereneBuilder
import FullereneWidget
import datetime
__addon_description__ = "Fullerene Builder"
__addon_version__ = "1.2"
__addon_date__ = datetime.date(2013, 10, 9)
__plugins__ = [FullereneBuilder.FullereneBuilder]

def reloadPlugins():
    reload(FullereneBuilder)
    reload(FullereneWidget)



