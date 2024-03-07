# Koester_DAmodels
white dwarf model spectra and an interpolator to easily generate white dwarf spectra


# Install
currently, you can only install this from github (or manually).

```
pip install git+https://github.com/janvanroestel/Koester_DAmodels
```

# Quickstart
You can import the module and make an instance of the WDmodels class.

```
import Koester_DAmodels.WDspec

DA = Koester_DAmodels.WDspec.WDmodels()

modelspectrum = DA.get_spectrum(T=12345.6,logg=7.89)
```

