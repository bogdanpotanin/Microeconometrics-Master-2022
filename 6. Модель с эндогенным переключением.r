# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 6. Модель с эндогенным переключением
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("mvtnorm")                       # симуляции из многомерного
                                         # нормального распределения

library("numDeriv")                      # численное дифференцирование

library("stringr")                       # работа со строками

library("tm")

library("switchSelection")               # метод Хекмана и модель с эндогенным
                                         # бинарным переключением

library("hpa")                           # моменты усеченного 
                                         # нормального распределения

#---------------------------------------------------
# Часть 1. Оценивание параметров
#---------------------------------------------------

set.seed(123)
# Зададим количество наблюдений
n <- 10000
# Истинные значения регрессионных коэффициентов
# каждого из уравнений
beta.1 <- c(1, 3, 6)
beta.0 <- c(-1, 2, 5)
gamma  <- c(0.5, 1, 2)
# Создадим независимые переменные из 
# многомерного нормального распределения
X <- rmvnorm(n,                                   # симулируем n наблюдений из многомерного
                                                  # нормального распределения
             c(0, 0, 0),                          # с нулевым вектором математических ожиданий и
             matrix(c(1, 0.2, 0.3,                # следующей ковариационной матрице
                      0.2, 1, -0.1,
                      0.3, -0.1, 1),
                    ncol = 3,
                    byrow = FALSE))
X[, 1] <- as.numeric(X[, 1] > 0.5)
# Симулируем случайные ошибик из двумерного
# нормального распределения
sigma.1 <- 3.5                                      # стандартное отклонение случайной
sigma.0 <- 3                                        # ошибки основного уравнения

rho.1 <- 0.8                                        # корреляция между случайными ошибками
rho.0 <- 0.5
rho.12 <- 0.6
epsilon_cov <- matrix(c(1,               rho.0 * sigma.0,            rho.1 * sigma.1,
                        rho.0 * sigma.0, sigma.0 ^ 2,                rho.12 * sigma.0 * sigma.1,
                        rho.1 * sigma.1, rho.12 * sigma.0 * sigma.1, sigma.1 ^ 2),
                      ncol = 3)
epsilon <- rmvnorm(n, c(0, 0, 0), epsilon_cov)    # случайные ошибки
# Латентная зависимая переменная
  # уравнения отбора 
z.star <- gamma[1] + gamma[2] * X[, 1] + 
          gamma[3] * X[, 2] + epsilon[, 1]
  # основного уравнения
y.star.0 <- beta.0[1] + beta.0[2] * X[, 2] + 
            beta.0[3] * X[, 3] + epsilon[, 2]
y.star.1 <- beta.1[1] + beta.1[2] * X[, 2] + 
            beta.1[3] * X[, 3] + epsilon[, 3]

# Наблюдаемая зависимая переменная
  # основного уравнения
z <- ifelse(z.star >= 0, 1, 0)
  # уравнения отбора 
y <- ifelse(z == 1, y.star.1, y.star.0)
# Посмотрим на результат
head(cbind(z, y, y.star.1, y.star.0))

# Сюжет:
# Представим, что вы изучаете, как уравнение
# зарплаты различается в зависимости от факта
# наличия у индивида высшего образования, при 
# это для простоты опустим неслучайный отбор
# в число занятых
h <- data.frame("wage" = y,                                # зарплата
                "educ" = z,                                # факт наличия высшего образования
                "parents" = X[, 1],                        # наличие высшего образования у родителей
                "age" = X[, 2],                            # возраст
                "health" = X[, 3])                         # состояние здоровья
head(h)

# Оценим параметры модели при помощи
  # МНК
model.ls.0 <- lm(wage ~ age + health,                      # раздельно оцениваем уравнения
                 data = h[h$educ == 0, ])                  # по выборкам из индивидов с
model.ls.1 <- lm(wage ~ age + health,                      # высшим образованием и без него
                 data = h[h$educ == 1, ])
coef.ls.0 <- coef(model.ls.0)                              # сохраняем коэффициенты раздельно
coef.ls.1 <- coef(model.ls.1)                              # для каждого из уравнений
  # Модель с эндогенным переключением 
  # на основе ММП
model.mle <- msel(                              
  formula   = educ ~ parents + age,                        # уравнение отбора
  formula2  = wage ~ age + health,                         # целевые уравнения
  data      = h,                                           # данные
  estimator = "ml",                                        # метод расчета ММП
  groups    = matrix(c(0, 1), ncol = 1),                   # возможные значения
                                                           # уравнения отбора
  groups2   = matrix(c(0, 1), ncol = 1))                   # возможные режимы                     
summary(model.mle)                                         # результат оценивания
coef.mle.0  <- coef(model.mle, type = "coef2", regime = 0) # сохраним оценки коэффициентов
coef.mle.1  <- coef(model.mle, type = "coef2", regime = 1)
cov.mle     <- coef(model.mle, type = "cov12")[[1]]        # оценки ковариаций между
cov.mle.0   <- cov.mle[1]                                  # случайными ошибками
cov.mle.1   <- cov.mle[2]                    
sigma.mle   <- sigma(model.mle)                            # стандартное отклонение
sigma.mle.0 <- sigma.mle[1]                                # случайной ошибки                          
sigma.mle.1 <- sigma.mle[2]                                
rho.mle.0   <- cov.mle.0 / sigma.mle.0
rho.mle.1   <- cov.mle.1 / sigma.mle.1
  # Модель с эндогенным переключением 
  # на основе двухшаговой процедуры
model.2st <- msel(                              
  formula   = educ ~ parents + age,                        # уравнение отбора
  formula2  = wage ~ age + health,                         # целевые уравнения
  data      = h,                                           # данные
  estimator = "2step",                                     # метод расчета ММП
  groups    = matrix(c(0, 1), ncol = 1),                   # возможные значения
                                                           # уравнения отбора
  groups2   = matrix(c(0, 1), ncol = 1))                   # возможные режимы                     
summary(model.2st)                                         # результат оценивания
coef.2st.0  <- coef(model.2st, type = "coef2", regime = 0) # сохраним оценки коэффициентов
coef.2st.1  <- coef(model.2st, type = "coef2", regime = 1)

# Сравним оценки и истинные значения
  # регрессионные коэффициенты
  # нулевого уравнения
data.frame("Real"          = beta.0,                       # истинные значения
           "Least Squares" = coef.ls.0,                    # МНК оценки
           "Switch MLE"    = coef.mle.0,                   # оценки ММП 
           "Switch 2step"  = coef.2st.0)                   # двухшаговые оценки
  # регрессионные коэффициенты
  # первого уравнения
data.frame("Real"          = beta.1,                       # истинные значения
           "Least Squares" = coef.ls.1,                    # МНК оценки
           "Switch MLE"    = coef.mle.1,                   # оценки ММП 
           "Switch 2step"  = coef.2st.1)                   # двухшаговые оценки
  # корреляция случайных ошибок
  # нулевого уравнения
data.frame("Real" = rho.0,                                 # истиннoе значение
           "Least Squares" = 0,                            # МНК оценки
           "Switch 2ST" = rho.2st.0,                       # оценки двухшаговые
           "Switch MLE" = rho.mle.0)                       # оценки ММП 
  # корреляция случайных ошибок
  # первого уравнения
data.frame("Real" = rho.1,                                 # истиннoе значение
           "Least Squares" = 0,                            # МНК оценки
           "Switch MLE" = rho.mle.1)                       # оценки ММП

#---------------------------------------------------
# Часть 2. Расчет предсказаний и предельных эффектов
#---------------------------------------------------

# Создадим индивида
Boris <- data.frame("parents" = 1,
                    "age"     = -0.3,
                    "health"  = 0.2)

# Рассчитаем оценку безусловного математического 
# ожидания зависимой переменной для каждого
# из уравнений, то есть E(y0*) и E(y1*)
  # E(y0*)
wage.star.0 <- predict(model.mle,                 # модель
                       newdata = Boris,           # данные
                       type    = "val",           # возвращаем матожидание
                       group2  = 0)               # режим 
  # E(y1*)
wage.star.1 <- predict(model.mle,                 # модель
                       newdata = Boris,           # данные
                       type    = "val",           # возвращаем матожидание
                       group2  = 1)               # режим 

# Рассчитаем оценку условного математического 
# ожидания зависимой переменной основного 
# уравнения, то есть E(y*|z)
  # E(y0|z = 0)
wage.0 <- predict(model.mle,                      # модель
                  newdata = Boris,                # данные
                  type    = "val",                # возвращаем матожидание
                  group2  = 1,                    # режим
                  group   = 0)                    # условие   
  # E(y1|z = 1)
wage.1 <- predict(model.mle,                      # модель
                  newdata = Boris,                # данные
                  type = "val",                   # возвращаем матожидание
                  group2 = 1,                     # режим
                  group = 1)                      # условие 

# Рассчитаем предельный эффект возраста
# на зарплату людей с высшим образованием
# то есть на E(y*|z = 1)
me_age <- predict(model.mle,                 # модель
                  newdata = Boris,           # данные
                  type    = "val",           # возвращаем матожидание
                  group2  = 1,               # режим
                  group   = 1,               # условие
                  me      = "age",           # переменная
                  test      = TRUE)          # стандартные ошибки
# Протестируем гипотезу о равенстве 
# предельного эффекта нулю
summary(me_age)

#---------------------------------------------------
# Часть 3. Пример оценивания
#---------------------------------------------------

# Загрузим данные
data(cps)
help(cps)
data <- cps

# Оценим уравнение зарплаты
model <- msel(formula  = work ~ age + nchild + health + bachelor + master,
              formula2 = lwage ~ age + bachelor + master,
              data     = data)
summary(model)

# Ослабим допущение о нормальном распределении
# случайной ошибки уравнения отбора
model2 <- msel(formula = work ~ age + nchild + health + bachelor + master,
               formula2 = lwage ~ age + bachelor + master,
               marginal = list(hpa = 3), data = data)
summary(model2)

# Проверим гипотезу о нормальности случайной ошибки
# уравннеия отбора
lrtest(model, model2)
c(parametric = AIC(model), semiparametric = AIC(model2))

# Ослабим допущение о совместном нормальном распределении
# и воспользуемся методом Ньюи
model3 <- msel(formula   =  work ~ age + nchild + health + bachelor + master,
               formula2  = lwage ~ age + bachelor + master,
               marginal  = list(logist = NULL),
               estimator = "2step", degrees = 2,
               data      = data)
summary(model3)

# Проверим гипотезу о том, что отдача от бакалаврского
# и магистерского уровня образования одинаковая
fn_test <- function(object)
{
  val <- coef(object, type = "coef2")[[1]][, "master"] -
         coef(object, type = "coef2")[[1]][, "bachelor"] 
  return (val)
}
test <- test_msel(model, fn = fn_test)
summary(test)

# Модель с эндогенным бинарным регрессором
model4 <- msel(formula  = list(work ~ age + nchild + health + 
                                      bachelor + master,
                               educ ~ age + sbachelor + smaster),
               formula2 = lwage ~ age + bachelor + master,
               data = data)
summary(model4)

# Модель зарплаты мужа и жены
model5 <- msel(formula  = list(work  ~ age + nchild + health + 
                                       bachelor + master,
                               swork ~ age + nchild + health + 
                                       sbachelor + smaster),
               formula2 = list(lwage  ~ age + bachelor + master,
                               slwage ~ age + sbachelor + smaster),
               data = data)
summary(model5)

# Модель с эндогенным переключением
groups  <- cbind(c(1, 1, 1, 0, 0, 0), c(0, 1, 2, 0, 1, 2))
groups2 <- matrix(c(0, 1, 2, -1, -1, -1), ncol = 1)
model6  <- msel(formula = list(work ~ age + nchild + health + 
                                      bachelor + master,
                               educ ~ age + sbachelor + smaster),
                formula2 = lwage ~ age + health,
                groups = groups, groups2 = groups2,
                data = data)
print(model6)
summary(model6)

# Посмотрим на структуру модели
struct_msel(model6)
