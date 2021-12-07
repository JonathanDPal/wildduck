fwhm_multiplier="$1"
python flux_measurement.py WildDuck R $fwhm_multiplier &
python flux_measurement.py WildDuck B $fwhm_multiplier &
python flux_measurement.py WildDuck V $fwhm_multiplier &
python flux_measurement.py m52 R $fwhm_multiplier &
python flux_measurement.py m52 B $fwhm_multiplier &
python flux_measurement.py m52 V $fwhm_multiplier &
