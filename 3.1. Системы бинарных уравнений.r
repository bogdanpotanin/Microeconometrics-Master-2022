# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 3.1. Системы бинарных уравнений
# --------

# Отключим scientific notation
options(scipen = 999)

library("switchSelection")                           # системы уравнений
library("mnorm")                                     # многомерное нормальное
                                                     # распределение
library("GJRM")                                      # системы уравнений с копулами

#---------------------------------------------------
# Симуляция данных
#---------------------------------------------------

# Для удобства представим, что мы симулируем процесс,
# определяющий дефолт по кредиту и факт наличия
# стабильной работы

# Симулируем данные
set.seed(123)                                              # для воспроизводимости
n          <- 10000                                        # число индивидов в выборке
h          <- data.frame(income = exp(rnorm(n, 10, 0.4)))  # доход
h$age      <- round(runif(n, 20, 100))                     # возраст
educ       <- t(rmultinom(n, 1, c(0.5, 0.3, 0.2)))         # уровень образования
h$educ_1   <- as.numeric(educ[, 1] == 1)                   # среднее образование
h$educ_2   <- as.numeric(educ[, 2] == 1)                   # среднее специальное образование
h$educ_3   <- as.numeric(educ[, 3] == 1)                   # высшее образование
h$credit   <- runif(n, 10000, 100000)                      # объем кредита
h$marriage <- rbinom(n, 1, 0.7)                            # официальный брак

# Симулируем случайные ошибки из двумерного
# нормального распределения
rho <- -0.5                                                # корреляция между
                                                           # случайными ошибками уравнений
eps <- rmnorm(n, 
              mean = c(0, 0), 
              sigma = matrix(c(1, rho, rho, 1), 
                             ncol = 2))
head(eps, 5)
eps1 <- eps[, 1]                                           # случайные ошибки уравнения дефолта
eps2 <- eps[, 2]                                           # случайные ошибки уравнения
                                                           # стабильной работы
     
# Зададим регрессоинные коэффициенты                                 
beta1 <- c( 3, -0.7, -0.02, 0.01,                        # оцениваемые регрессионные 
           -0.7, -0.9, 0.02,  -1)                          # коэффициенты уравнения дефолта
              
beta2 <- c(-2, 0.02, 2)                                    # оцениваемые регрессионные 
                                                           # коэффициенты уравнения
                                                           # стабильной работы

# Для уравнения стабильной работы создадим
  # линейный предиктор
stable_lp <- beta2[1] + beta2[2] * h$age + 
                        beta2[3] * h$marriage
  # латентную переменную
stable_star <- stable_lp + eps2 
  # зависимую переменную
h$stable <- as.numeric(stable_star >= 0)

# Для уравнения дефолта сформируем
  # линейный предиктор
default_lp <- beta1[1] + 
              beta1[2] * log(h$income) +
              beta1[3] * h$age +
              beta1[4] * (h$age ^ 2 / 100) +
              beta1[5] * h$educ_2 +
              beta1[6] * h$educ_3 +
              beta1[7] * sqrt(h$credit) +
              beta1[8] * h$stable
  # латентную переменную
default_star <- default_lp + eps1
  # зависимую переменную
h$default <- as.numeric(default_star >= 0)

# Посмотрим на маржинальные распределения
table(h$default)
table(h$stable)

# Посмотрим на совместное распределение
table(default = h$default, stable = h$stable)

# Итоговые данные
head(h, 10)

# Потенциальные исходы дефолта
default0 <- as.numeric((default_star - beta1[8] * h$stable) >= 0)
default1 <- as.numeric((default_star - beta1[8] * h$stable  + beta1[8]) >= 0)

# Средний эффект воздействия стабильной работы на дефолт
ATE0 <- mean(default1 - default0)

# Средний эффект воздействия стабильной работы на дефолт, 
# среди подвергшихся воздействию
ATET0 <- mean((default1 - default0)[h$stable == 1])

# --------------------------------------------
# Описание переменных:
# income    - доход
# age       - возраст
# educ_1    - среднее образование
# educ_2    - среднее специальное образование
# educ_3    - высшее образование
# credit    - объем кредита
# marriage  - состоит в браке
# default   - факт дефолта
# stable    - стабильная работа
# --------------------------------------------

#---------------------------------------------------
# Часть 1. Иерархическая система бинарных уравнений
#---------------------------------------------------

# Подготовим формулы
stable_formula  <- stable  ~ age + marriage
default_formula <- default ~ I(log(income)) + age + I(age ^ 2 / 100) +
                             educ_2 + educ_3 + 
                             I(sqrt(credit)) + stable

# Оценим параметры двумерной пробит модели
model <- msel(formula = list(stable_formula,
                             default_formula),
              data = h)                            
summary(model) 

# Достанем оценки коэффициентов
coef_default <- coef(model, type = "coef", eq = "default")
coef_stable  <- coef(model, type = "coef", eq = "stable")

# Попробуем оценить модель без учета эндогенности стабильной работы
# оценив обычную пробит модель (несостоятельная оценка)
model_ic        <- msel(formula = default_formula, data = h)      
coef_default_ic <- coef(model_ic, type = "coef", eq = "default")

# Сравним оценки и истинные значения
cbind(true = beta1[-1], bivprobit = coef_default, probit = coef_default_ic)
cbind(true = beta2[-1], estimate = coef_stable)

#---------------------------------------------------
# Часть 2. Оценивание вероятностей
#---------------------------------------------------

# Создадим отдельного индивида
Boris <- data.frame(age      = 35, income = 85000, 
                    educ_2   = 0,  educ_3 = 1,      
                    marriage = 1,  credit = 100000, 
                    stable   = 1)

# Оценим маржинальную вероятность дефолта Бориса 
# P(default = 1)
p_d1 <- predict(model, type = "prob", group = c(-1, 1), newdata = Boris)

# Оценим маржинальную отсутствия у Бориса стабильной работы 
# P(stable = 0)
p_s0 <- predict(model, type = "prob", group = c(0, -1), newdata = Boris)

# Оценим вероятность того, что у Бориса нет стабильной работы
# и произошел дефолт
# P(stable = 0, default = 1)
p_s0d1 <- predict(model, type = "prob", group = c(0, 1), newdata = Boris,
                  exogenous = list(stable = 0))

# Оценим вероятность того, что у Бориса произойдет дефолт, при условии,
# что у него нет стабильной работы
# P(default = 1 | stable = 0)
p_d1_s0 <- predict(model, type = "prob", group = c(0, 1), newdata = Boris,
                   given_ind = 1)

#---------------------------------------------------
# Часть 3. Оценивание предельных эффектов
#---------------------------------------------------

# Оценим предельный эффект возраста на P(default = 1)
me_p_d1 <- predict(model, type = "prob", group = c(-1, 1), newdata = Boris,
                   me = "age", test = TRUE)
summary(me_p_d1)

# Оценим предельный эффект брака на P(stable = 0)
me_p_s0 <- predict(model, type = "prob", group = c(0, -1), newdata = Boris,
                   me = "marriage", test = TRUE, eps = c(0, 1))
summary(me_p_s0)

# Оценим предельный эффект возраста на P(stable = 0, default = 1)
me_p_s0d1 <- predict(model, type = "prob", group = c(0, 1), newdata = Boris,
                     exogenous = list(stable = 0), me = "age", test = TRUE)
summary(me_p_s0d1)

# Оценим предельный эффект возраста на P(default = 1 | stable = 0)
me_p_d1_s0 <- predict(model, type = "prob", group = c(0, 1), newdata = Boris,
                      given_ind = 1, me = "age", test = TRUE,
                      exogenous = list(stable = 0))
summary(me_p_d1_s0)

#---------------------------------------------------
# Часть 4. Оценивание эффектов воздействия
#---------------------------------------------------

# Оценим средний предельный эффект стабильной работы на вероятность дефолта
# Важно: этот средний предельный эффект в данном случае совпадает со средним
#        эффектом воздействия стабильной работы на дефолт
ATE <- predict(model, type = "prob", group = c(-1, 1), 
               me = "stable", test = mean, eps = c(0, 1))
summary(ATE)

# Оценим средний эффект воздействия на вероятность дефолта, среди людей
# со стабильной работой
ATET_fn <- function(object)
{
  treated <-  object$data[object$data$stable == 1, ]
  p1 <- predict(object, type = "prob", group = c(1, 1), given_ind = 1,
                exogenous = list(stable = 1), newdata = treated)
  p0 <- predict(object, type = "prob", group = c(1, 1), given_ind = 1,
                exogenous = list(stable = 0), newdata = treated)
  ATET <- mean(p1 - p0)
  return(ATET)
}
ATET <- test_msel(model, ATET_fn)
summary(ATET)

# Альтернативная версия кода станет доступна
# после исправления бага в версии 2.0.1 пакета switchSelection
# ATET <- predict(model, type = "prob", group = c(1, 1), given_ind = 1,
#                 me = "stable", test = mean, eps = c(0, 1),
#                 newdata = model$data[model$data$stable == 1, ])
# summary(ATET)

# Сравним наши оценки со значениями, близкими к истинным
cbind(precise = ATE0, estimate  = ATE$tbl$val)
cbind(precise = ATET0, estimate = ATET$tbl$val)

#---------------------------------------------------
# Часть 5. Использование копул в системах
#          бинарных уравнений
#---------------------------------------------------

# Оценим еще две модели, предполагая, что у второй
# случайной ошибки логистическое распределение
# и рассматривая различные типы связи между случайными
# ошибками за счет копулы
  # Фрэнка
model_bp_1 <- gjrm(formula = list(default_formula,
                                  stable_formula),
                   data = h,
                   model = "B",
                   
                   margins = c("probit",                      # в первом уравнении распределение случайной
                               "logit"),                      # ошибки будет нормальным, а во втором - логистическим
                   copula = "F")                              # выбираем копулу Фрэнка
summary(model_bp_1)
  # Гумбеля
model_bp_2 <- gjrm(formula = list(default_formula,
                                  stable_formula),
                   data = h,
                   model = "B",
                   
                   margins = c("probit", "logit"),            # выбираем маржинальные распределения
                   copula = "G0")                             # выбираем копулу Гумбеля
summary(model_bp_2)

# Сравним модели по AIC
data.frame("Gaussian" = AIC(model),
           "Frank"    = AIC(model_bp_1),
           "Gumbel"   = AIC(model_bp_2))

# Таблицу с копулами и их описанием можно найти в:
# Marra G, Wyszynski K (2016), Semi-Parametric Copula Sample Selection Models for 
# Count Responses. Computational Statistics and Data Analysis, 104, 110-129.
# Другие работы авторов, связанные с пакетом GJRM
citation("GJRM")

#---------------------------------------------------
# Часть 5. Учет гетероскедастиночности и отклонений
#          от нормальности
#---------------------------------------------------

# Подключим данные
data("cps")
help(cps)
data <- cps

# Сформируем переменную на наличие высшего образования
# у женщины и супруга
data$higher  <- data$bachelor + data$master
data$shigher <- data$sbachelor + data$smaster

# Построим иерархическую модель связи между образование и работой
model.1 <- msel(list(higher ~ age + I(age ^ 2 / 100) + shigher,
                     work   ~ age + higher + nchild + health),
                data = data)
summary(model.1)

# Учтем гетероскедастиность
model.2 <- msel(list(higher ~ age + I(age ^ 2 / 100) + shigher | age,
                     work   ~ age + higher + nchild + health   | health),
                data = data)
summary(model.2)

# Рассмотрим полупараметрическую модель
model.3 <- msel(list(higher ~ age + I(age ^ 2 / 100) + shigher | age,
                     work   ~ age + higher + nchild + health   | health),
                marginal = list(logistic = NULL, hpa = 5), data = data)
summary(model.3)

# Сравним качество моделей
c(basic           = AIC(model.1), 
  heteroscedastic = AIC(model.2), 
  semiparametric  = AIC(model.3))

# Добавим регуляризацию для ковариации и коэффициента
# при высшем образовании в уравнении занятости
summary(model.1, show_ind = TRUE)
model.4 <- msel(list(higher ~ age + I(age ^ 2 / 100) + shigher,
                     work ~ age + higher + nchild + health),
                data = data, 
                regularization = list(ridge_ind = c(5, 8), 
                                      ridge_scale = 1000))
summary(model.4)

# Сформируем индивида
Alice <- data.frame(work   = 1, age = 30,   shigher = 0, 
                    higher = 1, nchild = 1, health  = 5)

# Посчитаем вероятность того, что Алиса работает и 
# имеет высшее образование
# P(higher = 1, work = 1)
p.11 <- predict(model.1, newdata = Alice, group = c(1, 1))

# Рассчитаем вероятность того, что Алиса работает
# P(work = 1)
p.x1 <- predict(model.1, newdata = Alice, group = c(-1, 1))

# Найдем вероятность того, что у Алисы есть высшее образование
# P(higher = 1)
p.1x <- predict(model.1, newdata = Alice, group = c(1, -1))

# Вычислим вероятность того, что Алиса работает, при условии, что
# у нее есть высшее образование (двумя способами)
# P(work = 1 | higher = 1)
pc <- predict(model.1, newdata = Alice, group = c(1, 1), given_ind = 1,
              exogenous = list(higher = 1))
pc <- p.11 / p.1x

# Оценим предельный эффект возраста на условную вероятность
# P(work = 1 | higher = 1)
me.age <- predict(model.1, newdata = Alice, group = c(1, 1), given_ind = 1,
                  exogenous = list(higher = 1), me = "age", test = TRUE)
summary(me.age)

# Оценим средний эффект воздействия высшего образования на занятость
ATE <- predict(model.1, type = "prob", group = c(-1, 1), 
               me = "higher", eps = c(0, 1), test = mean)
summary(ATE)
