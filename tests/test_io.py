import hytools as ht
from os.path import expanduser
from importlib import reload
import matplotlib.pyplot as plt
base = expanduser("~")

reload(ht)

# Open NEON test
neon = ht.openNEON("%s/data/neon.h5" % base)
neon.load_data()


red = neon.get_wave(660)
plt.matshow(red)


# Open AVIRIS test
aviris = ht.openENVI("%s/data/aviris" % base)
aviris.load_data()


red = aviris.get_wave(660)
plt.matshow(red)