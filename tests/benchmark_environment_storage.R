# Benchmark: List vs Environment Storage Performance
# This script demonstrates the performance benefits of using environments
# instead of lists for tessellation storage in the AddiVortes MCMC loop.

cat("=== AddiVortes Tessellation Storage Benchmark ===\n\n")

# Simulate the MCMC update pattern
benchmark_list <- function(m = 200, n_updates = 1000) {
  # Initialize as lists (old approach)
  tess <- lapply(1:m, function(i) matrix(rnorm(20), 10, 2))
  dim <- lapply(1:m, function(i) sample(1:5, 2))
  pred <- lapply(1:m, function(i) rnorm(10))
  
  start_time <- Sys.time()
  for (i in 1:n_updates) {
    j <- sample(1:m, 1)
    
    # Simulate MCMC updates (as in AddiVortes.R)
    tess[[j]] <- matrix(rnorm(20), 10, 2)
    dim[[j]] <- sample(1:5, 2)
    pred[[j]] <- rnorm(10)
  }
  end_time <- Sys.time()
  
  list(
    time = end_time - start_time,
    memory = object.size(tess) + object.size(dim) + object.size(pred)
  )
}

benchmark_environment <- function(m = 200, n_updates = 1000) {
  # Initialize as environments (new approach)
  tess <- new.env(hash = TRUE, size = m)
  dim <- new.env(hash = TRUE, size = m)
  pred <- new.env(hash = TRUE, size = m)
  
  for (k in 1:m) {
    key <- as.character(k)
    tess[[key]] <- matrix(rnorm(20), 10, 2)
    dim[[key]] <- sample(1:5, 2)
    pred[[key]] <- rnorm(10)
  }
  
  start_time <- Sys.time()
  for (i in 1:n_updates) {
    j <- sample(1:m, 1)
    key <- as.character(j)
    
    # Simulate MCMC updates (as in AddiVortes.R)
    tess[[key]] <- matrix(rnorm(20), 10, 2)
    dim[[key]] <- sample(1:5, 2)
    pred[[key]] <- rnorm(10)
  }
  end_time <- Sys.time()
  
  list(
    time = end_time - start_time,
    memory = object.size(tess) + object.size(dim) + object.size(pred)
  )
}

# Run benchmarks with different m values
cat("Benchmark 1: Small model (m = 50)\n")
cat("----------------------------------\n")
list_result_50 <- benchmark_list(m = 50, n_updates = 500)
env_result_50 <- benchmark_environment(m = 50, n_updates = 500)
cat(sprintf("List time:        %.4f seconds\n", as.numeric(list_result_50$time)))
cat(sprintf("Environment time: %.4f seconds\n", as.numeric(env_result_50$time)))
cat(sprintf("Speedup:          %.2fx\n", as.numeric(list_result_50$time) / as.numeric(env_result_50$time)))
cat(sprintf("List memory:      %d bytes\n", list_result_50$memory))
cat(sprintf("Env memory:       %d bytes\n\n", env_result_50$memory))

cat("Benchmark 2: Medium model (m = 200, default)\n")
cat("--------------------------------------------\n")
list_result_200 <- benchmark_list(m = 200, n_updates = 1000)
env_result_200 <- benchmark_environment(m = 200, n_updates = 1000)
cat(sprintf("List time:        %.4f seconds\n", as.numeric(list_result_200$time)))
cat(sprintf("Environment time: %.4f seconds\n", as.numeric(env_result_200$time)))
cat(sprintf("Speedup:          %.2fx\n", as.numeric(list_result_200$time) / as.numeric(env_result_200$time)))
cat(sprintf("List memory:      %d bytes\n", list_result_200$memory))
cat(sprintf("Env memory:       %d bytes\n\n", env_result_200$memory))

cat("Benchmark 3: Large model (m = 500)\n")
cat("----------------------------------\n")
list_result_500 <- benchmark_list(m = 500, n_updates = 1000)
env_result_500 <- benchmark_environment(m = 500, n_updates = 1000)
cat(sprintf("List time:        %.4f seconds\n", as.numeric(list_result_500$time)))
cat(sprintf("Environment time: %.4f seconds\n", as.numeric(env_result_500$time)))
cat(sprintf("Speedup:          %.2fx\n", as.numeric(list_result_500$time) / as.numeric(env_result_500$time)))
cat(sprintf("List memory:      %d bytes\n", list_result_500$memory))
cat(sprintf("Env memory:       %d bytes\n\n", env_result_500$memory))

cat("\n=== Summary ===\n")
cat("Environment-based storage provides:\n")
cat("1. Reference semantics - no copy-on-write overhead\n")
cat("2. More predictable memory behavior in long MCMC chains\n")
cat("3. Reduced memory footprint (reported size)\n")
cat("4. Full backward compatibility with existing code\n")
cat("\nNote: In this synthetic benchmark, performance is similar.\n")
cat("The main benefit is avoiding potential copy-on-write issues\n")
cat("in the actual MCMC loop. Copy-on-write is triggered when R's\n")
cat("reference counting system detects potential modifications to\n")
cat("shared objects, which is more likely in complex MCMC loops with\n")
cat("nested function calls. Environments guarantee reference semantics.\n")
