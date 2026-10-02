program test_volcano
  use patmo, only: patmo_loadInitialProfile
  use patmo_photo, only: loadPhotoMetric
  use patmo_volc, only: volcano_prerun_settings, patmo_volc_configurePreRun, patmo_volc_runPreRun
  implicit none
  type(volcano_prerun_settings)::config

  call patmo_volc_configurePreRun(config)
  call loadPhotoMetric(trim(config%photoMetricFile))
  call patmo_loadInitialProfile(trim(config%profileFile),unitH="km",unitX="1/cm3")
  call patmo_volc_runPreRun(config)
end program test_volcano
