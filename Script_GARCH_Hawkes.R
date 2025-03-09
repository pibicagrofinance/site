

# Loading R packages

library(quantmod)
library(dplyr)
library(ggplot2)
library(tidyr)
library(timetk)
library(patchwork)
library(plotly)
library(rugarch)

########################################
# 1) Carregar e Preparar os Dados
########################################

# Baixar dados históricos do milho (ZC=F) do Yahoo Finance
getSymbols("ZC=F", src = "yahoo", from = "2020-01-01", to = Sys.Date())

# Criar data.frame com a data e o preço de fechamento
milho_data <- data.frame(
  Date  = index(`ZC=F`),
  Close = as.numeric(Cl(`ZC=F`))
)

# Calcular os log-retornos diários e remover linhas com NA
milho_data <- milho_data |>
  mutate(Return = c(NA, diff(log(Close)))) |>
  na.omit()

########################################
# 2) Estimar GARCH(1,1) com Distribuição Skewed Student-t
########################################

# Vetor de retornos
ret <- milho_data$Return

# Especificação do modelo GARCH(1,1) com distribuição "sstd" (skewed Student-t)
spec <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "sstd"  # skewed Student-t
)

# Ajuste do modelo
fit <- ugarchfit(
  spec   = spec,
  data   = ret,
  solver = "hybrid"
)

# Exibe o sumário do modelo
show(fit)

# Extrai a volatilidade condicional estimada
cond_vol <- sigma(fit)

# Converte a volatilidade para vetor numérico e anexa ao data.frame original
milho_data$Volatility <- as.numeric(cond_vol)

########################################
# 3) Reestruturar os Dados para Plot com Facetas
########################################

# Seleciona as séries de interesse e converte para o formato "long"
milho_data_long <- milho_data |>
  select(Date, Close, Return, Volatility) |>
  pivot_longer(
    cols = c(Close, Return, Volatility),
    names_to = "Serie",
    values_to = "Valor"
  )

# Define a ordem e os labels para as facetas: Preços, Log-Retornos e Volatilidade
milho_data_long$Serie <- factor(milho_data_long$Serie,
                                levels = c("Close", "Return", "Volatility"),
                                labels = c("Preços", "Log-Retornos", "Volatilidade"))

########################################
# 4) Plotar as Séries Empilhadas com timetk
########################################

plot_faceted <- milho_data_long |>
  group_by(Serie) |>
  plot_time_series(
    .date_var    = Date,
    .value       = Valor,
    .interactive = FALSE,   # Defina TRUE se preferir interatividade
    .facet_ncol  = 1,       # Uma série por linha (empilhadas verticalmente)
    .smooth      = FALSE,
    .title       = "Preços, Log-Retornos e Volatilidade Condicional do Milho"
  ) +
  theme(strip.background = element_rect(fill = "white", colour = "white"))

# Exibe o gráfico interativo
ggplotly(plot_faceted)

########################################
# 5) Calcular Resíduos Padronizados e Definir Eventos Extremos
########################################

# Calcular os resíduos padronizados: ret_std = Return / Volatility
milho_data <- milho_data |>
  mutate(Standardized = Return / Volatility)

# Definir threshold para eventos extremos (ex.: |resíduo| > 2)
threshold_std <- 2
extreme_idx <- which(abs(milho_data$Standardized) > threshold_std)
# Os retornos extremos (magnitude dos jumps) serão usados para os saltos
extreme_returns <- milho_data$Return[extreme_idx]

# Converter as datas dos eventos extremos para tempo (em dias) relativo ao primeiro dia
first_date <- min(milho_data$Date)
event_times_hist <- as.numeric(milho_data$Date[extreme_idx] - first_date)

# Tempo total de observação (em dias)
T_obs <- as.numeric(max(milho_data$Date) - first_date)

########################################
# 6) Calibrar o Processo de Hawkes com Eventos Extremos Históricos
########################################

# Função de log-verossimilhança negativa para um Hawkes com kernel exponencial
neg_log_lik <- function(params, event_times, T_total) {
  mu_h    <- params[1]
  alpha_h <- params[2]
  beta_h  <- params[3]
  
  # Penaliza parâmetros não positivos
  if(mu_h <= 0 || alpha_h <= 0 || beta_h <= 0) return(1e10)
  
  n <- length(event_times)
  log_sum <- 0
  for(i in seq_along(event_times)) {
    ti <- event_times[i]
    if(i == 1) {
      sum_exp <- 0
    } else {
      sum_exp <- sum(exp(-beta_h * (ti - event_times[1:(i-1)])))
    }
    lambda_ti <- mu_h + alpha_h * sum_exp
    log_sum <- log_sum + log(lambda_ti)
  }
  integral_term <- mu_h * T_total + (alpha_h / beta_h) * sum(1 - exp(-beta_h * (T_total - event_times)))
  ll <- log_sum - integral_term
  return(-ll)  # retorna o negativo para minimização
}

# Chute inicial para os parâmetros do Hawkes
init_par <- c(mu_h = 0.1, alpha_h = 0.5, beta_h = 1.0)

# Otimização para calibrar o modelo Hawkes
res_hawkes <- optim(
  par         = init_par,
  fn          = neg_log_lik,
  event_times = event_times_hist,
  T_total     = T_obs,
  method      = "L-BFGS-B",
  lower       = c(1e-6, 1e-6, 1e-6),
  upper       = c(Inf, Inf, Inf)
)

# Parâmetros calibrados do Hawkes
mu_h_hat    <- res_hawkes$par[1]
alpha_h_hat <- res_hawkes$par[2]
beta_h_hat  <- res_hawkes$par[3]

cat("Parâmetros Hawkes Estimados:\n")
cat("mu =", mu_h_hat, "\n")
cat("alpha =", alpha_h_hat, "\n")
cat("beta =", beta_h_hat, "\n")

########################################
# 7) Simular 3 Trajetórias dos Preços para os Próximos 252 Dias
########################################

# Utilizar a previsão do GARCH para obter a média e volatilidade basal futura
horizon <- 252
garch_forecast <- ugarchforecast(fit, n.ahead = horizon)
mu_forecast    <- as.numeric(fitted(garch_forecast))    # média prevista (geralmente próxima de zero)
sigma_forecast <- as.numeric(sigma(garch_forecast))       # volatilidade prevista

# Preço inicial: último preço histórico
P0 <- tail(milho_data$Close, 1)

# Função para simular eventos futuros via Hawkes (algoritmo de thinning)
simulateHawkes <- function(mu_h, alpha_h, beta_h, T_start, T_end, history = NULL) {
  if(is.null(history)) history <- numeric(0)
  current_events <- history
  t <- T_start
  new_events <- c()
  
  while(t < T_end) {
    # Cota superior para a intensidade baseada em eventos recentes (últimos 10 dias)
    n_recent <- sum(current_events > (t - 10))
    lambda_bar <- mu_h + alpha_h * n_recent
    w <- rexp(1, rate = lambda_bar)
    t_candidate <- t + w
    if(t_candidate > T_end) break
    sum_exp <- if(length(current_events[current_events < t_candidate]) > 0)
      sum(exp(-beta_h * (t_candidate - current_events[current_events < t_candidate]))) else 0
    lambda_tc <- mu_h + alpha_h * sum_exp
    if(runif(1) <= lambda_tc / lambda_bar) {
      new_events <- c(new_events, t_candidate)
      current_events <- c(current_events, t_candidate)
    }
    t <- t_candidate
  }
  return(new_events)
}

n_sim <- 3  # Número de trajetórias
simulated_paths <- list()

set.seed(123)  # Para reprodutibilidade

for(sim in 1:n_sim) {
  # Simular eventos extremos futuros via Hawkes para o intervalo [T_obs, T_obs+horizon]
  new_event_times <- simulateHawkes(mu_h_hat, alpha_h_hat, beta_h_hat,
                                    T_start = T_obs,
                                    T_end   = T_obs + horizon,
                                    history = event_times_hist)
  
  # Converter tempos simulados para dias relativos (1 a horizon)
  event_days <- floor(new_event_times - T_obs) + 1
  
  # Para cada dia, se houver eventos, somar os jumps.
  # Os valores dos jumps são amostrados dos retornos extremos históricos.
  jumps <- rep(0, horizon)
  if(length(event_days) > 0) {
    for(day in unique(event_days)) {
      n_events <- sum(event_days == day)
      jump_vals <- sample(extreme_returns, n_events, replace = TRUE)
      jumps[day] <- sum(jump_vals)
    }
  }
  
  # Simular os retornos basais para cada dia (usando a previsão do GARCH)
  baseline_returns <- rnorm(horizon, mean = mu_forecast, sd = sigma_forecast)
  
  # O retorno total é a soma do retorno basal e dos jumps
  total_returns <- baseline_returns + jumps
  
  # Simular a trajetória dos preços: P_t = P0 * exp(cumsum(total_returns))
  price_path <- P0 * exp(cumsum(total_returns))
  
  simulated_paths[[sim]] <- data.frame(
    Day = 1:horizon,
    Price = price_path,
    Return = total_returns,
    Trajectory = paste0("Sim", sim)
  )
}

# Combinar as trajetórias para plotagem
sim_df <- do.call(rbind, simulated_paths)

# Adicionar uma coluna 'Date' para as simulações (forecast)
last_date <- tail(milho_data$Date, 1)
sim_df <- sim_df |>
  mutate(Date = as.Date(last_date) + Day)

########################################
# 8) Combinar Histórico e Forecast no Mesmo Gráfico
########################################

# Data frame com preços históricos (usando a série 'Close')
historical_df <- milho_data |>
  select(Date, Close) |>
  rename(Price = Close)

# Plot: Linha sólida para dados históricos, linhas tracejadas para simulações,
# e uma linha vertical indicando o início do forecast.
p_sim <- ggplot() +
  geom_line(data = historical_df, aes(x = Date, y = Price),
            color = "black", size = 1) +
  geom_line(data = sim_df, aes(x = Date, y = Price, color = Trajectory),
            linetype = "3313") +
  geom_vline(xintercept = as.numeric(last_date),
             linetype = "dotted", color = "red") +
  labs(title = "Preços do Milho: Histórico e Simulação de 3 Trajetórias (Próximos 252 Dias)",
       x = "Data", y = "Preço") +
  theme_minimal()

ggplotly(p_sim)




