# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 6.2. Модель Тобина
# --------

# Для удобства отключим 
# экспоненциальную запись чисел
options(scipen = 999)

# Подключим дополнительные библиотеки
library("crch")                          # регрессия с усечением

library("hpa")                           # моменты усеченного 
                                         # нормального распределения

# Симулируем данные
library("mvtnorm")

# Симулируем данные
set.seed(777)
n <- 2000
h <- data.frame(inc     = exp(rnorm(n, 10, 0.7)),     # доход
                weather = pmin(rpois(n, 0.5) + 1, 8), # погодные условия (1 - плохие, 8 - лучшие)
                deliv   = round(runif(n, 0, 127)))    # цена доставки

# Случайная ошибка
sigma   <- sqrt(700000)
epsilon <- rnorm(n, mean = 0, sd = sigma)

# Не усеченная зависимая переменная
beta       <- c(3500, 0.07, -1000, -30)
spend_star <- beta[1]             + 
              beta[2] * h$inc     +
              beta[3] * h$weather +
              beta[4] * h$deliv   + 
              epsilon

# Точки усечения
tr_left <- 2000

# Создаем усечение
h$spend                     <- spend_star
h$spend[h$spend <= tr_left] <- tr_left
nrow(h)
head(h, 5)

#---------------------------------------------------
# Часть 1. Оценивание параметров модели Тобина
#---------------------------------------------------

# Оценим параметры модели Тобина
tr <- 2000                                           # точка усечения
model_tb <- crch(spend ~ inc + weather + deliv,      # формула
                 data = h,                           # данные                       
                 left = tr,                          # нижнее (левое) усечение
                 truncated = FALSE,                  # модель Тобина
                 dist = "gaussian")                  # распределение случайной
                                                     # ошибки
summary(model_tb)                                    # посмотрим результат

# Сохраним оценки коэффициентов
est_tb   <- coef(model_tb)                           # достанем оценки
coef_tb  <- est_tb[-length(est_tb)]                  # оценки регрессионных
                                                     # коэффициентов
sigma_tb <- exp(est_tb[length(est_tb)])              # достаем оценку 
                                                     # стандартного отклонения

# Воспользуемся усеченной регрессией для сравнения   # точка усечения
model_tr <- crch(spend ~ inc + weather + deliv,      # формула
                 data = h[h$spend > tr_left, ],      # данные                       
                 left = tr,                          # нижнее (левое) усечение
                 truncated = TRUE,                   # усеченная регрессия 
                 dist = "gaussian")                  # распределение случайной
                                                     # ошибки
summary(model_tr)                                    # посмотрим результат

# Сохраним оценки коэффициентов
est_tr   <- coef(model_tr)                           # достанем оценки
coef_tr  <- est_tr[-length(est_tr)]                  # оценки регрессионных
                                                     # коэффициентов
sigma_tr <- exp(est_tr[length(est_tr)])              # достаем оценку 
                                                     # стандартного отклонения

# Воспользуемся обычным МНК для сравнения
model_lm <- lm(spend ~ inc + weather + deliv,
               data = h)
coef_lm  <- coef(model_lm)
sigma_lm <- sigma(model_lm)

# Сравним истинные коэффициенты с оценками
cbind(true  = beta, 
      tobin = coef_tb,
      trunc = coef_tr, 
      ls    = coef_lm)

#---------------------------------------------------
# Часть 2. Расчет вероятностей и предельных эффектов
#---------------------------------------------------

# Рассчитаем оценку безусловного математического
# расходов для конкретного индивида
Boris <- data.frame(inc     = 25000,
                    weather = 3,
                    deliv   = 1)

# Предскажем ожидаемые расходы без 
# учета цензурирования, то есть E(spend*)
spend_est <- predict(model_tb, newdata = Boris)

# Оценим вероятность того, что Борис
# расходует больше 2000
prob_est <- 1 - pnorm((tr - spend_est) / sigma_tb)

# Предскажем условное математическое
# ожидание пасходов
# вычислим E(e | spend* >= tr_left)
  # с помощью функции
epsilon_E <- truncatedNormalMoment(k = 1,                      # момент
                                   x_lower = tr - spend_est,   # нижнее усечение
                                   x_upper = Inf,              # верхние усечение
                                   mean = 0,                   # математическое ожидание
                                   sd = sigma_tb)              # стандартное отклонение
  # вручную
a <- (tr - spend_est) / sigma_tb
lambda <- dnorm(a) / (1 - pnorm(a))
epsilon_E <- sigma_tb * lambda
  # посчитаем E(spend* | spend* >= tr_left)
spend_est_cond <- spend_est + epsilon_E
  # рассчитаем E(spend) = P(spend* >= tr_left) * E(spend* | spend* >= tr_left)
  #                       P(spend* < tr_left) * E(spend | spend* < tr_left) =
  #                       P(spend* >= tr_left) * E(spend* | spend* >= tr_left) +
  #                       (1 - P(spend* < tr_left)) * tr_left
spend_est_cens <- prob_est * spend_est_cond + (1 - prob_est) * tr_left

# Рассмотрим предельные эффекты на
  # E(spend* | spend* >= tr_left)
ME_cond <- coef_tb * (1 - lambda * (lambda - a))
  # E(spend)
ME_cens <- coef_tb * prob_est

#---------------------------------------------------
# Часть 3. Гетероскедастичность в модели Тобина
#---------------------------------------------------

# Допустим, что дисперсия случайной ошибки
# зависит от цены доставки и погоды
set.seed(321)
gamma   <- c(log(sigma), -0.1, 0.01)                       # коэффициенты уравнения
                                                           # дисперсии
sigmai  <- exp(gamma[1] +                                  # различныедисперсии
               h$weather * gamma[2] + 
               h$deliv * gamma[3])
epsiloni <- sigmai * rnorm(n)                              # гетероскедостичные случайные
                                                           # ошибки
h$spend <- (spend_star - epsilon) + epsiloni               # учтем гетероскедастичность
h$spend[h$spend < tr_left] <- tr_left

model_htobit <- crch(spend ~ inc + weather + deliv |       # основное уравнение
                                   weather + deliv,        # уравнение дисперсии
                     data       = h,                       # данные
                     left       = tr_left,                 # усечение
                     link.scale = "log")                   # тип уравнения дисперсии           
summary(model_htobit) 

# Аргумент dist позволяет оценивать параметры
# модели при альтернативном допущении
# о распределении случайных ошибок
model_htobit_t <- crch(spend ~ inc + weather + deliv |
                                     weather + deliv,      
                       data       = h,                        
                       left       = tr_left,                  
                       link.scale = "log",      
                       dist       = "student")             # будем использовать
                                                           # распределение Стьюдента
summary(model_htobit_t) 

# Встроенные функции для расчета вероятностей 
# работает по аналогии

# Некоторые дополнительные полезные функции
logLik(model_tb)                                           # правдоподобие
AIC(model_tb)                                              # AIC
BIC(model_tb)                                              # BIC

# Сравним модели
c(normal = AIC(model_htobit), student = AIC(model_htobit_t))

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 3.1.    Сравните модели с нормальным распределением
#         случайных ошибок и распределением Стьюдента
#         по критерию AIC
# 3.2*.   Проверьте, что произойдет, если в уравнение
#         дисперсии будет добавлена незначимая 
#         переменная
# 3.3**.  Для индивида с произвольными характеристиками
#         оцените предельный эффект возраста на:
#         1)     E(y*)
#         2)     P(y > upper)
#         3)     P(y < lower)
#         4)     E(y|lower < y < upper)
#         5)     E(y)
#         6*)    E(y|y > 0)

#---------------------------------------------------
# Часть 4. Самостоятельная реализация модели Тобина
#---------------------------------------------------

# Уберем гетероскедастичность
h$spend <- spend_star
h$spend[h$spend < tr_left] <- tr_left

# Запишем функцию правдоподобия
# модели Тобина
lnL <- function(par, y, X, 
                tr_left = -Inf, 
                tr_right = Inf)
{
  par_n <- length(par)                                     # число оцениваемых
                                                           # параметров
  sigma <- par[par_n]              
  beta <- matrix(par[-par_n], ncol = 1)
  
  if (sigma <= 0)                                          # следим, чтобы в оптимизатор
  {                                                        # не было подано отрицательное
    return(-(1e+100))                                      # стандартное отклонение
  }
  
  X <- cbind(1, X)                                         # добавляем константу
  
  XB <- X %*% beta                                         # считаем произведение
                                                           # регрессоров на коэффициенты
  n <- nrow(X)                                             # число наблюдений
  
  L <- rep(NA, n)                                          # вклады в правдоподобие
  
  cond_1 <- (y == tr_left)                                 # различные условия
  cond_2 <- (y == tr_right)                                # дают различные вклады
  cond_3 <- ((y > tr_left) & (y < tr_right))               # в функцию правдоподобия

  if (!is.infinite(tr_left))
  {
    L[cond_1] <- pnorm(tr_left - XB[cond_1, ],             # считаем вклады в
                       sd = sigma)                         # функцию правдоподобия
  }
  
  if (!is.infinite(tr_right))
  {
    L[cond_2] <- pnorm(XB[cond_2, ] - tr_right, 
                       sd = sigma)
  }
  
  L[cond_3] <- dnorm(y[cond_3] - XB[cond_3, ], 
                     sd = sigma)
  
  lnL.val <- sum(log(L))                                   # считаем логарифм
                                                           # функции правдоподобия
  
  return(lnL.val)
}
# Осуществим оптимизацию
x0 <- c(coef_tr, sigma_tr)                                 # в качестве начальной
                                                           # точки возьмем оценки
                                                           # усеченной регрессии
opt_mle <- optim(par = x0, fn = lnL,                       # запускаем оптимизатор
                 X = cbind(h$inc, h$weather, h$deliv), 
                 y = h$spend,
                 tr_left = tr_left,
                 method = "Nelder-Mead", 
                 control = list(maxit = 100000000,          
                                fnscale = -1,
                                reltol = 1e-10,
                                abstol = 1e-10))
x1 <- opt_mle$par                                          # сохраняем оценки
# Сравним истинные и прогнозные значения
data.frame("beta" = beta, 
           "beta_hat" = x1[-length(x1)])

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 4.1*.   Перепишите функцию правдоподобия, используя
#         репараметризацию Олсена, а именно, замены:
#         beta* = beta / sigma
#         sigma* = 1 / sigma
# 4.2**.  Оцените гипотезы о значимости параметров
#         используя оценки, полученные с использованием
#         репараметризации Олсена
# 4.3*.   Реализуйте модель Тобина со случайными ошибками,
#         имеющими распределение Стьюдента с 5-ю степенями
#         свободы