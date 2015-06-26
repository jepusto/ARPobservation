PND <- function(data, phase, base_phase, increase = TRUE) {
  if (!increase) data <- -1 * data
  100 * mean(data[phase != base_phase] > max(data[phase == base_phase]))
}


PEM <- function(data, phase, base_phase, increase = TRUE) {
  if (!increase) data <- -1 * data
  med <- median(data[phase == base_phase])
  100 * mean((data[phase != base_phase] > med) + 0.5 * (data[phase != base_phase] == med))
}


PAND <- function(data, phase, base_phase, increase = TRUE) {
  if (!increase) data <- -1 * data
  m <- sum(phase == base_phase)
  n <- sum(phase != base_phase)
  X <- sort(data[phase == base_phase])
  Y <- sort(data[phase != base_phase])
  ij <- expand.grid(i = 1:m, j = 1:n)
  ij$no_overlap <- mapply(function(i, j) X[i] < Y[j], i = ij$i, j = ij$j)
  ij$overlap <- with(ij, i + n - j + 1)
  overlaps <- with(ij, max(overlap * no_overlap))
  
  100 * overlaps / length(data)
}


IRD <- function(data, phase, base_phase, increase = TRUE) {
  pand <- PAND(data = data, phase = phase, base_phase = base_phase, increase = increase) / 100
  m <- sum(phase == base_phase)
  n <- sum(phase != base_phase)
  
  ((m+n)^2 * pand - m^2 - n^2) / (2 * m * n)
}


NAP <- function(data, phase, base_phase, increase = TRUE) {
  if (!increase) data <- -1 * data
  XY <- expand.grid(x = data[phase==base_phase], y = data[phase!=base_phase])
  100 * mean(with(XY, (y > x) + 0.5 * (y == x)))
}


Tau <- function(data, phase, base_phase, increase = TRUE) {
  nap <- NAP(data = data, phase = phase, base_phase = base_phase, increase = increase)
  2 * nap / 100 - 1
}

SMD <- function(data, phase, base_phase, ...) {
  treat_phase <- levels(phase)[(base_phase != levels(phase))]
  y_bar <- tapply(data, phase, mean)[c(base_phase, treat_phase)]
  s_sq <- tapply(data, phase, var)[c(base_phase, treat_phase)]
  n <- table(phase)[c(base_phase, treat_phase)]
  s_pooled <- sqrt(sum((n - 1) * s_sq) / sum(n - 1))
  (y_bar[2] - y_bar[1]) / s_pooled
}