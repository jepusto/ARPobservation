context("Testing the AIR observation filter")

test_that("AIR data is consistent with component procedures", {
  BS <- r_behavior_stream(n = 100, mu = 3, lambda = 10, 
                          F_event = F_exp(), F_interim = F_exp(), stream_length = 1000)
  interval_length <- 10
  rest_length <- 5
  AIR <- augmented_recording(BS, interval_length, rest_length)
  MTS <- momentary_time_recording(BS, interval_length, summarize = FALSE)
  PIR <- interval_recording(BS, interval_length, rest_length, partial = TRUE, summarize = FALSE)
  WIR <- interval_recording(BS, interval_length, rest_length, partial = FALSE, summarize = FALSE)
  
  expect_true(all(AIR[,"MTS",] == MTS))
  expect_true(all(AIR[-1,"PIR",] == PIR))
  expect_true(all(AIR[-1,"WIR",] == WIR))
})
