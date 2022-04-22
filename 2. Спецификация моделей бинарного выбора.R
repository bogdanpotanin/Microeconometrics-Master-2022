# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 2. Спецификация моделей бинарного выбора
# --------

# Отключим scientific notation
options(scipen = 999)

# Симулируем данные, содержащую информацию
# о характеристиках заемщика, а также о том,
# наступил ли у него дефолт по ипотеке.
# Все переменные измерены в условных единицах
set.seed(12345)                                            # для воспроизводимости
n <- 10000                                                 # число индивидов в выборке
h <- data.frame(ind = rep(1:n))                            # датафрейм для хранения данных
h$inc <- runif(n, 0, 1)                                    # доход в условных единицах
h$pay <- runif(n, 0, 1)                                    # ежемесячный платеж
h$age <- runif(n, 0, 1)                                    # возраст
h$ins <- rbinom(n, 1, 0.7)                                 # факт наличия страховки
h$chl <- rbinom(n, 1, 0.6)                                 # факт наличия детей
eps <- rnorm(n)                                            # случайная ошибка 

beta <- c(0.6, -3, 2, -5, 3.5, -0.8, 0, 0.8)               # оцениваемые регрессионные 
                                                           # коэффициенты


def_li <- beta[1] +                                        # линейный индекс,
  beta[2] * h$inc +                                        # отражающий вклад наблюдаемых
  beta[3] * h$pay +                                        # факторов в вероятность дефолта
  beta[4] * h$age +
  beta[5] * h$age ^ 2 +                                            
  beta[6] * h$ins +                          
  beta[7] * h$chl +
  beta[8] * (h$chl * h$inc)

def_star <- def_li + eps                                   # латентная переменная,
                                                           # отражающая склонность
                                                           # к дефолту
h$def <- as.numeric(def_star >= 0)                         # наблюдаемое значение переменной
mean(h$def)                                                # доля дефолтов

h$ind <- NULL                                              # уберем индексы
head(h, 5)                                                 # посмотрим на данные

#---------------------------------------------------
# Часть 1. Учет гетероскедастичности в моделях
#          бинарного выбора
#---------------------------------------------------

# Кратко о предполагаемом процессе генерации данных
# в пробит модели с гетероскедастичной случайной
# ошибкой:
# y_latent_i = x_i * beta + e_i          основное уравнение
# e_i ~ N(0, sigma_i ^ 2)                случайная ошибка
# sigma_i = h(w_i * tau)                 уравнение дисперсии
# w_i * tau                              линейный индекс
#                                        уравнения дисперсии
# w_i                                    независимые переменные
#                                        влияющие на дисперсию
# h(0) = 1                               обычно накладываемые
# h'(0) != 0                             условия на функцию h()
# h(t) = exp(t)                          некоторые примеры
#        |t + 1|                         функции h(), удовтелвтяорющей
#        (t + 1) ^ 2                     данным условиям
# В данной модели в качестве оцениваемых параметров
# выступают векторы коэффициентов beta и tau

# Внесем некоторые корректировки в процесс 
# генерации данных

# Симулируем гетероскедастичные случайные ошибки, 
# предполагая, что их дисперсия может зависеть от
# платежей и факта наличия детей
set.seed(123)
tau <- c(0.6, -0.3)
eps.var <- exp(tau[1] * h$pay +                           # дисперсия случайной ошибки
               tau[2] * h$chl) ^ 2                        # зависит от некоторых регрессоров
head(eps.var, 5)
eps <- rnorm(n, mean = 0, sd = sqrt(eps.var))             # случайные ошибки с
                                                          # различными дисперсия

                                                          # Симулируем зависимую переменную
  def_star <- def_li + eps                                  # латентная переменная,
                                                          # отражающая склонность
                                                          # к дефолту
h$def <- as.numeric(def_star >= 0)                        # наблюдаемое значение переменной
mean(h$def)                                               # доля дефолтов

# Оценим параметры пробит модели с 
# гетероскедастичной случайной ошибкой

library("glmx")                                           # пакет, позволяющий оценивать пробит
                                                          # модель с гетероскдестичной 
                                                          # случайной ошибкой
                                                          # Оценим пробит модель без
                                                          # учета гетероскедастичности
library("lmtest")                                         # тестирование гипотез
model.probit <- glm(formula = def ~ inc + pay +
                                    age + I(age ^ 2) +
                                    ins + chl +
                                    chl * inc,
                    data = h,                                  
                    family = binomial(link = "probit")) 
summary(model.probit)

# Оценим пробит модель без
# с учетом гетероскедастичности
model.hetprobit <- hetglm(formula = def ~ inc + pay +      # линейный индекс
                                          age +            # основного уравнения
                                          I(age ^ 2) +
                                          ins + chl +
                                          chl * inc |
                                          pay + chl,       # линейный индекс
                                                           # уравнения дисперсии
                          data = h,                                 
                          family = binomial(link = "probit"),
                          link.scale = "log")
summary(model.hetprobit)
# В функции hetglm() link.scale указывает, 
# в каком виде представлена ошибка
# Имеются следующие варианты:
# 1. identity  ---  sigma_i        =  w_i * tau
# 2. log       ---  log(sigma_i)   =  w_i * tau  =>  sigma_i = exp(w_i * tau_i)
# 3. sqrt      ---  sqrt(sigma_i)  =  w_i * tau  =>  sigma_i = (w_i * tau_i) ^ 2

# Достанем полученные оценки
beta.est <- model.hetprobit$coefficients$mean              # оценки коэффициентов при переменных
# основного уравнения
tau.est <- model.hetprobit$coefficients$scale              # оценки коэффициентов при переменных
# в уравнении дисперсии

# Сравним истинные значения и оценки
# коэффициентов линейного индекса
# уравнения дисперсии
rbind(true = tau,
      est = tau.est)

# Достанем оценки стандартных отклонений
# случайных ошибок
sigma.est <- predict(model.hetprobit, type = "scale")
head(sigma.est, 5)

# Осуществим тест на гомоскедастичность:
# H0: tau = 0
lrtest(model.hetprobit, model.probit)

# Предскажем вероятности и 
# скорректированный линейный индекс
prob.est <- predict(model.hetprobit, type = "response")    # вероятности
head(prob.est, 5)
y.li.est <- predict(model.hetprobit, type = "link")        # скорректированный 
head(y.li.est, 5)                                          # линейный индекс

# Рассчитаем предельные эффекты
# для индивида
Boris <- data.frame(inc = 0.2,                             # характеристики Бориса
                    pay = 0.1,
                    age = 0.3,
                    ins = 1,
                    chl = 1)

# Предварительные расчеты
prob.Boris <- predict(model.hetprobit, newdata = Boris,    # оценка вероятности 
                      type = "response")                   # дефолта Бориса
li.Boris.adj <- predict(model.hetprobit, newdata = Boris,  # оценка отношения линейного
                        type = "link")                     # индекса Бориса к стнадртному
                                                           # отклонению случайно ошибки
sigma.Boris <- predict(model.hetprobit, newdata = Boris,   # оценка стандартного
                       type = "scale")                     # отклонения случайной
                                                           # ошибки Бориса
li.Boris <- li.Boris.adj * sigma.Boris                     # оценка линейного
                                                           # индекса Бориса

# Используем встроенную функцию
library("margins")
ME.Boris <- margins(model.hetprobit, 
                    data = Boris)
summary(ME.Boris)
ME.pay <- ME.Boris$dydx_pay

# Считаем предельный эффект аналитически   
ME.pay.1 <- dnorm(li.Boris, sd = sigma.Boris) * 
  (beta.est["pay"] - 
   li.Boris * tau.est["pay"])

# Считаем предельный эффект с помощью
# численного дифференцирования
delta <- 1e-6                                              # приращение                                                 
Boris.delta <- Boris
Boris.delta$pay <- Boris$pay + delta                       # приращение по возрасту
prob.Boris.delta <- predict(model.hetprobit, 
                            newdata = Boris.delta,
                            type = "response")
ME.pay.2 <- (prob.Boris.delta - prob.Boris) / delta

# Сравним посчитанные предельные эффекты
rbind(margins = ME.pay,
      analytical = ME.pay.1,
      numeric = ME.pay.2)

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 1.1. Используя встроенные данные Mroz87 из библиотеки
#      sampleSelection и пробит модель с гетероскедастичной
#      случайно ошибкой определите, как на вероятность
#      занятости (lfp) влияют возраст (age), образование (educ),
#      факт проживания в городе (city) и число несовершеннолетних
#      детей (kids5 и kids618). При этом предполагается, что
#      гетероскедастичность может быть обусловлена возрастом
#      и уровнем образования. Далее, для 28-летнего индивида 
#      без высшего образования и с доходом 20000 оцените:
#      1)    вероятность занятости
#      2)    предельный эффект возраста на вероятность занятости
#      3)    предельный эффект проживания в городе на вероятность занятости
#      4*)   предельный эффект возраста на вероятность занятости, если
#            возраст входит в линейный индекс квадратично
#      5)    повторите предыдущие пункты, используя различные подходы
#            к определению формы уравнения дисперсии: см. аргумент link.scale
#      6**)  стандартную ошибку оценки вероятности занятости

#---------------------------------------------------
# Часть 2. Проверка гипотезы о гомоскедастичности
#          при помощи LM теста
#---------------------------------------------------

library("numDeriv")                                        # библиотека для
                                                           # численного дифференцирования
                                                           # Оценим модель с гетероскедастичной
                                                           # случайной ошибкой
model.hetprobit <- hetglm(formula = def ~ inc + pay +
                                          age + 
                                          I(age ^ 2) +
                                          ins + chl +
                                          I(chl * inc) |
                                          pay + chl,
                          data = h,                                 
                          family = binomial(link = "probit"))
beta.est <- model.hetprobit$coefficients$mean
tau.est <- model.hetprobit$coefficients$scale

# Проверим гипотезу о гомоскедастичности
# при помощи LM теста, преимущество которого
# заключается в том, что нет необходимости
# предполагать конкретную форму функции h(t):
# достаточно наложить на нее лишь пару ограничений:
# 1) h(0) = 1                              
# 2) h'(0) != 0                             
HetprobitLnL <- function(x,                                # коэффициенты
                         y,                                # зависима переменная
                         X,                                # регрессоры основного уравнения
                         W,                                # регрессоры уравнения дисперсии
                         scale_fn = exp,                   # функция уравнения дисперсии h()
                         is_aggregate = TRUE)              # возвращаем функцию правдоподобия (TRUE)   
  # или отдельные вклады (FALSE)          
{
  m_X <- ncol(X)
  m_W <- ncol(W)
  
  beta <- matrix(x[1:m_X], ncol = 1)                       # вектор beta коэффициентов и
  tau <- matrix(x[(m_X + 1):(m_X + m_W)], ncol = 1)        # вектор дополнительных параметров  
                                                           # переводим в матрицу с одним столбцом
  
  y_li_mean <- X %*% beta                                  # оценка линейного индекса
  y_li_scale <- W %*% tau                                  # латентной переменной
  y_li_scale_fn <- scale_fn(y_li_scale)
  
  n_obs <- nrow(X)                                         # количество наблюдений
  
  L_vec <- matrix(NA, nrow = n_obs,                        # вектор столбец вкладов наблюдений
                  ncol = 1)                                # в функцию правдоподобия
  
  is_y_0 <- (y == 0)                                       # вектор условий y = 0
  is_y_1 <- (y == 1)                                       # вектор условий y = 1
  
  L_vec[is_y_1] <- pnorm(y_li_mean[is_y_1], 
                         sd = y_li_scale_fn[is_y_1])       # вклад наблюдений для которых yi = 1
  L_vec[is_y_0] <- 1 - pnorm(y_li_mean[is_y_0],
                             sd = y_li_scale_fn[is_y_0])   # вклад наблюдений для которых yi = 0
  
  lnL_vec <- log(L_vec)                                    # логарифмы вкладов
  
  if(!is_aggregate)                                        # возвращаем вклады
  {                                                        # при необходимости
    return(lnL_vec)
  }
  
  lnL <- sum(lnL_vec)                                      # логарифм функции правдоподобия
  
  return(lnL)
}

# Достанем данные
df.hetprobit <- model.frame(model.hetprobit)               # все регрессоры
X.mat <- cbind(1, as.matrix(df.hetprobit[                  # регрессоры основного
  names(df.hetprobit) %in% names(beta.est)]))              # уравнения
head(X.mat , 5)
W.mat <- as.matrix(df.hetprobit[                           # регрессоры уравнения
  names(df.hetprobit) %in% names(tau.est)])                # дисперсии
head(W.mat , 5)
y.vec <- df.hetprobit[, 1]                                 # зависимая переменная

# Достанем оценки ограниченной модели
x.est.R <- c(model.probit$coefficients, 
             rep(0, ncol(W.mat)))
n.R <- length(x.est.R)                                     # добавим имена
names(x.est.R)[(n.R - 1):n.R] <- paste(colnames(W.mat),    # для красоты
                                       "sigma")
print(x.est.R)

# Рассчитаем правдоподобие полной модели в точке,
# определяемой оценками, полученными по ограниченной
# модели
lnL.R <- HetprobitLnL(x.est.R, 
                      y = y.vec,
                      X.mat, W.mat, exp)
lnL.R.grad <- grad(func = HetprobitLnL,                    # считаем градиент данной функции
                   x = x.est.R,                            # численным методом
                   y = y.vec, 
                   X = X.mat, W = W.mat,
                   scale_fn = exp)                         # замените exp на function(x)
                                                           #                 {
                                                           #                   return(abs(x + 1)})
                                                           #                 }
                                                           # и убедитесь, что результат не изменится
lnL.R.grad <- matrix(lnL.R.grad, ncol = 1)                 # градиент как матрица с одним столбцом
lnL.R.Jac <- jacobian(func = HetprobitLnL,                 # оцениваем асимптотическую ковариационную
                      x = x.est.R,                         # матрицу при помощи Якобиана, расcчитанного
                      y = y.vec,                           # численным методом, поскольку численно
                      X = X.mat, W = W.mat,                # рассчитать Гессиан достаточно точным 
                      scale_fn = exp,                      # образом не получается
                      is_aggregate = FALSE,
                      method.args = list(r = 8))

# Реализуем тест
LM.value <- t(lnL.R.grad) %*%                              # считаем статистику теста
  solve(t(lnL.R.Jac) %*% lnL.R.Jac) %*%                    # множителей Лагранжа
  lnL.R.grad
p.value <- 1 - pchisq(LM.value, df = ncol(W.mat))          # рассчитываем p-value теста

# Аналитический подход
li.R <- X.mat %*% x.est.R[1:ncol(X.mat)]  
jac.gamma <- (dnorm(li.R) / 
              pnorm((2 * y.vec - 1) * li.R)) * 
             (1 - 2 * y.vec) * li.R
jac.gamma <- as.vector(jac.gamma) * W.mat
jac.beta <- (dnorm(li.R) / pnorm((2 * y.vec - 1) * 
             li.R)) * (2 * y.vec - 1)
jac.beta <- as.vector(jac.beta) * X.mat
jac <- cbind(jac.beta, jac.gamma)
gr <- colSums(jac)
LM.value.2 <- gr %*% solve(t(jac) %*% jac) %*% gr
p.value.2 <- 1 - pchisq(LM.value.2, df = ncol(W.mat))

#---------------------------------------------------
# Часть 3. Тестирование гипотезы о нормальном
#          распределении случайных ошибок
#---------------------------------------------------

library("numDeriv")                                        # численно дифференцирование

# Оценим пробит модель
model.probit <- glm(formula = def ~ inc + pay +            # оцениваем модель
                                    age + I(age ^ 2) +
                                    ins + chl +
                                    I(chl * inc),
                    data = h,                                  
                    family = binomial(link = "probit")) 

# Запишем функцию правдоподобия
# для модели со случайно ошибкой
# из распределения Пирсона
ProbitLnLExtended <- function(par,                         # вектор значений параметров
                              y,                           # зависимая переменная 
                              X,                           # матрица независимых переменных
                              is_aggregate = TRUE)         # при TRUE возвращаем логарифм
                                                           # функции правдоподобия, а при
                                                           # FALSE возвращаем вектор вкладов
{
  beta <- matrix(par[-c(1, 2)], ncol = 1)                  # вектор beta коэффициентов и
  theta <- matrix(par[c(1, 2)], ncol = 1)                  # вектор дополнительных параметров  
                                                           # переводим в матрицу с одним столбцом
  y_li <- X %*% beta                                       # оценка линейного индекса
  y_est <- y_li + theta[1] * y_li ^ 2 +                    # оценка математического ожидания 
           theta[2] * y_li ^ 3                             # латентной переменной
  
  n_obs <- nrow(X)                                         # количество наблюдений
  
  L_vec <- matrix(NA, nrow = n_obs,                        # вектор столбец вкладов наблюдений
                  ncol = 1)                                # в функцию правдоподобия
  
  is_y_0 <- (y == 0)                                       # вектор условий (y = 0)
  is_y_1 <- (y == 1)                                       # вектор условий (y = 1)
  
  L_vec[is_y_1] <- pnorm(y_est[is_y_1])                    # вклад наблюдений для которых yi = 1
  L_vec[is_y_0] <- 1 - pnorm(y_est[is_y_0])                # вклад наблюдений для которых yi = 0
  
  lnL_vec <- log(L_vec)                                    # логарифмы вкладов
  
  if(!is_aggregate)                                        # возвращаем вклады
  {                                                        # при необходимости
    return(lnL_vec)
  }
  
  lnL <- sum(lnL_vec)                                      # логарифм функции правдоподобия
  
  return(lnL)
}
# Воспользуемся созданной функцией
# Оценки модели при справедливом ограничении,
# накладываемом нулевой гипотезой
beta.est <- coef(model.probit)                             # достаем оценки из обычной пробит
beta.R <- c(0, 0, beta.est)                                # модели и приравниваем значения
names(beta.R)[c(1, 2)] <- c("theta1", "theta2")            # дополнительных параметров к значениям,
                                                           # предполагаемым нулевой гипотезой
print(beta.R)
# Создадим матрицу регрессоров
X.mat <- as.matrix(model.frame(model.probit))              # достаем датафрейм с регрессорами и
X.mat[, 1] <- 1                                            # первращаем его в матрицу, а также
colnames(X.mat)[1] <- "Intercept"                          # заменяем зависимую переменную на константу
head(X.mat, 5)
# Применим функцию
lnL.R <- ProbitLnLExtended(beta.R, h$def, X.mat)           # считаем логарифм функции правоподобия
                                                           # при ограничениях, совпадающую с логарифмом
# функции правдоподобия обычной пробит модели
lnL.R.grad <- grad(func = ProbitLnLExtended,               # считаем градиент данной функции
                   x = beta.R,                             # численным методом
                   y = h$def, 
                   X = X.mat)
lnL.R.grad <- matrix(lnL.R.grad, ncol = 1)                 # градиент как матрица с одним столбцом
lnL.R.Jac <- jacobian(func = ProbitLnLExtended,            # считаем Якобин данной функции
                      x = beta.R,                          # численным методом
                      y = h$def, 
                      X = X.mat,
                      is_aggregate = FALSE)
# Реализуем тест
LM.value.1 <- t(lnL.R.grad) %*%                            # считаем статистику теста
  solve(t(lnL.R.Jac) %*% lnL.R.Jac) %*%                    # множителей Лагранжа
  lnL.R.grad
p.value_1 <- 1 - pchisq(LM.value.1, df = 2)                # рассчитываем p-value теста
                                                           # множителей Лагранжа

# С использованием регрессии на единицы

# Достанем датафрейм, содержащий
# переменные модели
d <- model.frame(model.probit)                             # все переменные

# Рассчитаем предварительные величины
y.li.est <- predict(model.probit)                          # оценка линейного индекса                                     
F.est <- pnorm(y.li.est)                                   # функции от линейного           
f.est <- dnorm(y.li.est)                                   # индекса

# Вычислим обобщенные остатки
gr <- ((d[, 1] - F.est) /                                  # обобщенный остаток
         (F.est * (1 - F.est))) * f.est

# Считаем производные по коэффициентам
d_beta <- apply(X.mat, 2, function(x)                      # производные по
{                                                          # регресcионным коэффициентам
  x * gr
})
d_t1 <- (gr * y.li.est ^ 2)                                # производная по t1
d_t2 <- (gr * y.li.est ^ 3)                                # производная по t2

# Сравним аналитические и численные производные
grad_df <- data.frame("Numeric" = lnL.R.grad,
                      "Analytical" = colSums(cbind(d_t1, 
                                                   d_t2, 
                                                   d_beta)))
rownames(grad_df) <- c("t1", "t2", colnames(X.mat))
print(grad_df)

# Проводим LM тест
n <- nrow(d)                                               # число наблюдений
LM_df <- data.frame("my_ones" = rep(1, n),                 # вектор из единиц 
                    "d_" = d_beta,
                    d_t1, d_t2)          
head(LM_df, 5)
ones_regression <- summary(lm(my_ones ~. + 0,              # регрессия на вектор единиц без константы
                              data = LM_df))       
R2 <- ones_regression$r.squared                            # коэффициент детерминации регрессии
LM_value_2 <- R2 * n                                       # LM статистика
p.value_2 <- 1 - pchisq(q = LM_value_2, 
                        df = 2)

# Сравним полученные результаты и убедимся,
# что они полностью совпадают
c(LM.value.1, LM_value_2)                                  # сравниваем статистики
c(p.value_1, p.value_2)                                    # сравниваем p-value

#---------------------------------------------------
# Часть 4. Метод Галланта и Нички
#---------------------------------------------------

library("hpa")                                             # полу-непараметрический подход
library("ggplot2")                                         # красивые графики
library("Ecdat")                                           # встроенные данные

# По мотивам статьи:
# Parametric and Semi-Parametric Estimation of the 
# Binary Response Model of Labour Market

# Рассмотрим пример на данных, отражающих 
# занятость индивидов
help(Participation)

# Загрузим данные
data("Participation")
h <- Participation
h$lfp <- as.numeric(h$lfp == "yes")
h$foreign <- as.numeric(h$foreign == "yes")
h$age <- h$age * 10

# Пробит модель, описывающая занятость
model_pr <- glm(formula = lfp ~ lnnlinc +                  # логарифм нетрудового дохода
                                age + I(age ^ 2) +         # возраст   
                                educ +                     # образование в годах
                                nyc +                      # к-во маленьких детей
                                noc +                      # к-во взрослых детей
                                foreign,                   # иностранец 
                data = h,                                     
                family = binomial(link = "probit"))  
summary(model_pr)
coef_pr <- coef(model_pr)

# Как и в оригинальной работе будем использовать
# полином третьей степени
model_hpa <- hpaBinary(formula = lfp ~ I(-lnnlinc) +       # нормализуем коэффициент при 
                                       age + I(age ^ 2) +  # регрессоре lnnlinc к -1
                                       educ +                      
                                       nyc +                       
                                       noc +                       
                                       foreign,                   
                       data = h, 
                       K = 3,                              # степень полинома
                       cov_type = "sandwich")              # тип ковариационной матрицы
summary(model_hpa)
coef_hpa <- model_hpa$coefficients                         # достанем оценки коэффициентов

# Визуализируем результат
plot(model_hpa)                                            # оценка функции плотности
                                                           # случайной ошибки
# Сравним модели по AIC
AIC(model_hpa)    
AIC(model_pr)

# сравним оценки со стандартизированными 
# коэффициентами пробит модели
data.frame("Galland.Nychka" = coef_hpa, 
           "Probit" = coef_pr[-1] / (-coef_pr[2]))

# Рассмотрим пример на данных, отражающих готовность
# людей платить за сохранение парка
help(Kakadu)

# Загрузим данные
data("Kakadu")
h <- Kakadu
h$wtp <- as.numeric(h$answer != "nn")                      # переменная принимает значение 1,
                                                           # если индивид готов заплатить за
                                                           # сохранение парка больше некоторой суммы

# Модель, описывающая готовность индивида
# заплатить более некоторой суммы
model_pr <- glm(formula = wtp ~ mineparks +                # открытие производства в парковых
                                                           # зонах существенно уменьшает
                                                           # привлекательность парка
                                age +                      # возраст
                                sex +                      # пол (мужчина)
                                income +                   # доход в тысячах долларов
                                moreparks +                # нужно больше парков
                                wildlife +                 # важно сохранять дикую природу
                                aboriginal +               # важно учитывать интересы 
                                                           # коренных жителей
                                finben,                    # важно руководствоваться соображениями
                                                           # финансовой выгоды при использовании
                                                           # природных ресурсов
                data = h,                                    
                family = binomial(link = "probit"))          
summary(model_pr)

# Применим полунепараметрическую модель
model_hpa <- hpaBinary(formula = formula(model_pr),                   
                       data = h, 
                       K = 6,                              # степень полинома
                       cov_type = "sandwich")              # тип ковариационной матрицы
summary(model_hpa)
# ВАЖНО:
# В зависимости от типа ковариационной матрицы
# можно получить различный подход к интерпретации
# оценок:
# cov_type = "sandwich" - квази-максимальное правдоподобие
# cov_type = "bootstrap" - полу-непараметрический подход
# cov_type = "hessian" - параметрический подход


# Визуализируем результат
plot(model_hpa)

# Сравним модели по AIC
AIC(model_pr)
AIC(model_hpa)

# Сравним предсказанные вероятности
p_pr <- predict(model_pr, type = "response")
p_hpa <- predict(model_hpa, is_prob = TRUE)
head(cbind(probit = p_pr, GN = p_hpa))
plot(p_pr, p_hpa,                                          # визуально сравнение оценок
     xlab = "Gallant and Nychka", ylab = "Probit")         # вероятностей

# Сравним предельные эффекты пробит модели
# и полученные при помощи метода Галланта и Нички
ME_hpa <- model_hpa$marginal_effects                       # Галланта и Ничка
ME_probit <- margins(model_pr, type = "response")          # Пробит модель
plot(ME_hpa[, "age"], ME_probit$dydx_age,                  # визуально сравнение оценок
     xlab = "Gallant and Nychka", ylab = "Probit")         # предельного эффекта возраста

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 4.1.    Подберите оптимальную степень полинома K
#         руководствуясь соображениями минимизации
#         информационного критерия Акаике
# 4.2.    Сравните оценки асимптотической ковариационной
#         матрицы, получаемые различными методами
# 4.3.    Используя функцию dhpa аналитически оцените:
#         1*)   вероятность того, что индивид с вашими
#               характеристиками будет готов платить
#               за сохранение парка
#         2**)  предельный эффект возраста на
#               соответствующую вероятность