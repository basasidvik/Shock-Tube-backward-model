Shock-tube initial-state ML  -  order in which to run the files
================================================================

Put ALL .m files in one folder and make it MATLAB's current folder.
Delete old .mat / .csv files from earlier runs first (including sensor_model.mat).
All settings (cases, times, sensors, noise, floors, model) live in
shock_tube_config.m - edit that file, not the step scripts.

RUN, IN THIS ORDER (press Run on each):

  step1_generate_datasets.m     generate + store 10 FULL .mat files
                                (t = 0.04 ... 0.40, 50,000 cases, 256 points each)
  step2_export_sensor_csv.m     use only t = 0.2 and 0.4; write labelled CSVs
                                of the 6 sensor positions
  check_stored_data.m           (optional) confirm what is stored where
  step3_train_model.m           train the neural network -> sensor_model.mat
  step4_refine_with_physics.m   sharpen the NN estimate with the exact solver

  run_shock_tube.m              (optional, any time) solve + plot one case

Re-run only step2 (then 3 and 4) if you change the sensor positions, the sensor
noise or the times used; step1 only needs repeating if the cases, grid or the
stored times change.

Called automatically - do NOT press Run on these (they need inputs):
  shock_tube_config, generate_dataset_multitime, export_sensor_csv,
  riemann_exact, riemann_sample_batch, riemann_star, pressure_function,
  predict_initial_state, refine_initial_state
