# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 7. Модель с частотными данными
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("MASS")                                            # отрицательная биномиальная
                                                           # регрессия
library("pscl")                                            # модель с инфляцией нулей

# Будем моделировать посещение индивидом врача
set.seed(123)
n <- 10000
income <- runif(n)
age <- runif(n)
beta <- c(0.1, 0.8, 1)
lambda <- exp(beta[1] + beta[2] * income + beta[3] * age)

# Пуассоновская регрессия
doctor <- rpois(n, lambda)
model.pois <- glm(doctor ~ income + age, family = "poisson")
summary(model.pois)
# предсказания
lambda_est <- exp(predict(model.pois))
head(cbind(lambda, lambda_est))
# предельные эффекты
me_age <- lambda_est * coef(model.pois)["age"]
ame_age <- mean(me_age)
me_income <- lambda_est * coef(model.pois)["income"]
ame_income <- mean(me_income)

# Отрицательная биномиальная регрессия
theta <- 0.3
eps <- log(rgamma(n, shape = 1 / theta, scale = theta))
lambda <- exp(beta[1] + beta[2] * income + beta[3] * age + eps)
doctor <- rpois(n, lambda)
model.nb <- glm.nb(doctor ~ income + age)
summary(model.nb)
# предсказания
lambda_est <- exp(predict(model.nb))
head(cbind(lambda, lambda_est))
# Вместо theta пакет оценивает 1 / theta
theta_hat <- 1 / model.nb$theta
# Сравнение с моделью Пуассона
model.pois <- glm(doctor ~ income + age, family = "poisson")
lambda_est_pois <- exp(predict(model.pois))
head(cbind(lambda, lambda_est, lambda_est_pois))
# Проверка нулевой гипотезы о том, что
# можно пользоваться моделью Пуассона (нет overdispersion)
sd_theta_hat <- model.nb$SE.theta * abs(1 / model.nb$theta ^ 2)
z <- theta_hat / sd_theta_hat
p_value <- 2 * min(pnorm(z), 1 - pnorm(z))
# предельные эффекты
me_age <- lambda_est * coef(model.nb)["age"]
ame_age <- mean(me_age)
me_income <- lambda_est * coef(model.nb)["income"]
ame_income <- mean(me_income)

# Модель с инфляцией нулей
city <- rbinom(n, 1, 0.5)
gamma <- c(0.1, 3, -2, 1)
u <- rnorm(n)
z.li <- gamma[1] + gamma[2] * income + gamma[3] * age + gamma[4] * city
z.star <- z.li + u
z <- rep(1, n)
z[z.star < 0] <- 0
doctor[z == 0] <- 0
data <- data.frame(doctor, income, age, city)
# но в пакете оценивает, меняя местами 0 и 1, поэтому
# в уравнении нулей коэффициенты будут с обратными знаками
model.zim <- zeroinfl(doctor ~ income + age | 
                               income + age + city, 
                      data = data,
                      dist = "negbin", link = "probit")
summary(model.zim)
# вместо теты вновь оценивается 1 / theta
theta_hat <- 1 / model.zim$theta
# предсказания
lambda_est <- predict(model.zim)
head(cbind(lambda, lambda_est, z), 10)
# предельные эффекты
gamma_est <- -model.zim$coefficients$zero
beta_est <- model.zim$coefficients$count
li <- gamma_est["age"] * age + gamma_est["income"] * income +
      gamma_est["city"] * city
me_age <- lambda_est * pnorm(li) * beta_est["age"] +
          lambda_est * dnorm(li) * gamma_est["age"]
ame_age <- mean(me_age)

# Ручная реализация Пуассоновской регрессии
lnL <- function(par, X, y)
{
  beta <- par
  X <- cbind(1, X)
  li <- X %*% beta
  lambda <- exp(li)
  log_lik <- sum(dpois(x = y, lambda = lambda, log = TRUE))
  return(log_lik)
}
X <- cbind(income, age)
n_par <- ncol(X) + 1
y <- doctor
opt <- optim(par = rep(0, 3),                   # начальная точка
             fn = lnL,                          # максимизируемая функция правдоподобия
             method = "Nelder-Mead",            # оптимизационный алгоритм
             hessian = TRUE,                    # возвращаем Гессиан, он нам понадобится
             control = list(maxit = 10000,      # максимальное количество итераций
                            fnscale = -1,       # ставим -1 чтобы сделать задачу максимизационной
                            reltol = 1e-10),    # условие остановки алгоритма (termination condition)
             X = X, y = y)
beta_est <- opt$par
model.pois <- glm(doctor ~ income + age, family = "poisson")
beta_est_1 <- coef(model.pois)
rbind(manual = beta_est, automatic = beta_est_1)

# Модифицированная Пуассоновская регрессия
# пример с плохой идентифицируемостью
  # генерация данных
set.seed(123)
n <- 10000
income <- runif(n)
age <- runif(n)
beta <- c(0.1, 0.8, 1)
li <- beta[1] + beta[2] * income + beta[3] * age
tau <- 0.3
lambda <- abs(li) + tau * exp(li)
doctor <- rpois(n, lambda)
  # оценивание
lnL <- function(par, X, y)
{
  X <- cbind(1, X)
  beta <- par[1:ncol(X)]
  tau <- par[ncol(X) + 1]
  li <- X %*% beta
  lambda <- abs(li) + tau * exp(li)
  log_lik <- sum(dpois(x = y, lambda = lambda, log = TRUE))
  return(log_lik)
}
X <- cbind(income, age)
n_par <- ncol(X) + 2
y <- doctor
opt <- optim(par = rep(0.01, n_par),            # начальная точка
             fn = lnL,                          # максимизируемая функция правдоподобия
             method = "Nelder-Mead",            # оптимизационный алгоритм
             hessian = TRUE,                    # возвращаем Гессиан, он нам понадобится
             control = list(maxit = 10000,      # максимальное количество итераций
                            fnscale = -1,       # ставим -1 чтобы сделать задачу максимизационной
                            reltol = 1e-10),    # условие остановки алгоритма (termination condition)
             X = X, y = y)
beta_est <- opt$par[1:(n_par - 1)]
tau_est <- opt$par[n_par]
rbind(beta, beta_est)
c(tau, tau_est)
