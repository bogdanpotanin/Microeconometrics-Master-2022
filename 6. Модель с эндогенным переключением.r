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

library("sampleSelection")               # метод Хекмана

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
beta.1 <- c(1, 3, 5)
beta.0 <- c(-1, 2, 6)
gamma <- c(0.5, 1, 2)
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
sigma.1 <- 3                                      # стандартное отклонение случайной
sigma.0 <- 3.5                                    # ошибки основного уравнения

rho.1 <- 0.8                                      # корреляция между случайными ошибками
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
head(cbind(z, y))

# Сюжет:
# Представим, что вы изучаете, как уравнение
# зарплаты различается в зависимости от факта
# наличия у индивида высшего образования, при 
# это для простоты опустим неслучайный отбор
# в число тех, кто имеет высшее образование
h <- data.frame("wage" = y,                                # зарплата
                "educ" = z,                                # факт наличия высшего образования
                "parents" = X[, 1],                        # наличие высшего образования у родителей
                "age" = X[, 2],                            # возраст
                "health" = X[, 3])                         # состояние здоровья
head(h)

# Оценим параметры модели при помощи
  # МНК
model.ls.0 <- lm(wage ~ age + health,                      # раздельно оцениваем уравнения
               data = h[h$educ == 0, ])                    # по выборкам из индивидов с
model.ls.1 <- lm(wage ~ age + health,                      # высшим образованием и без него
                 data = h[h$educ == 1, ])
coef.ls.0 <- coef(model.ls.0)                              # сохраняем коэффициенты раздельно
coef.ls.1 <- coef(model.ls.1)                              # для каждого из уравнений
  # Модель с эндогенным переключением 
  # на основе ММП
model.mle <- selection(                              
  selection = educ ~ parents + age,                        # уравнение отбора
  outcome = list(wage ~ age + health,                      # целевые уравнения
                 wage ~ age + health),
  data = h,                                                # данные
  method = "ml",                                           # метод расчета ММП
  type = "5")                                              # укзываем модель с эндогенным переклюечением 
summary(model.mle)                                         # результат оценивания
coef.mle <- coef(model.mle, part = "outcome",              # сохраним оценки коэффициентов
                  prefix = TRUE)        
coef.mle.0 <- coef.mle[1:3]
coef.mle.1 <- coef.mle[4:6]
rho.mle.0 <- model.mle$estimate["rho1"]                    # оценка корреляции между
rho.mle.1 <- model.mle$estimate["rho2"]                    # случайными ошибками

sigma.mle.0 <- model.mle$estimate["sigma1"]                # стандартное отклонение
sigma.mle.1 <- model.mle$estimate["sigma2"]                # случайной ошибки

# Модель с эндогенным переключением 
# на основе двухшаговой процедуры
model.2st <- selection(                              
  selection = educ ~ parents + age,                        # уравнение отбора
  outcome = list(wage ~ age + health,                      # целевые уравнения
                 wage ~ age + health),
  data = h,                                                # данные
  method = "2step",                                        # метод расчета ММП
  type = "5")                                              # укзываем модель с эндогенным переклюечением         
summary(model.2st)                                         # результат оценивания
coef.2st <- coef(model.2st, part = "outcome",              # сохраним оценки коэффициентов
                 prefix = TRUE)        
coef.2st.0 <- coef.2st[1:3]
coef.2st.1 <- coef.2st[5:7]

# Сравним оценки и истинные значения
  # регрессионные коэффициенты
  # нулевого уравнения
data.frame("Real" = beta.0,                                # истинные значения
           "Least Squares" = coef.ls.0,                    # МНК оценки
           "Heckman MLE" = coef.mle.0,                     # оценки ММП 
           "Heckman 2step" = coef.2st.0)                   # двухшаговые оценки
  # регрессионные коэффициенты
  # первого уравнения
data.frame("Real" = beta.1,                                # истинные значения
           "Least Squares" = coef.ls.1,                    # МНК оценки
           "Heckman MLE" = coef.mle.1,                     # оценки ММП 
           "Heckman 2step" = coef.2st.1)                   # двухшаговые оценки
  # корреляция случайных ошибок
  # нулевого уравнения
data.frame("Real" = rho.0,                                 # истиннoе значение
           "Least Squares" = 0,                            # МНК оценки
           "Heckman MLE" = rho.mle.0)                      # оценки ММП 
  # корреляция случайных ошибок
  # первого уравнения
data.frame("Real" = rho.1,                                 # истиннoе значение
           "Least Squares" = 0,                            # МНК оценки
           "Heckman MLE" = rho.mle.1)                      # оценки ММП 

#---------------------------------------------------
# Часть 2. Расчет предсказаний и предельных эффектов
#---------------------------------------------------

# Создадим индивида
Boris <- data.frame("wage" = 1,
                    "educ" = 1,
                    "parents" = 1,
                    "age" = -0.3,
                    "health" = 0.2)

# Рассчитаем оценку безусловного математического 
# ожидания зависимой переменной для каждого
# из уравнений, то есть E(y0*) и E(y1*)
wage.star <- predict(model.mle, 
                     newdata = Boris, 
                     part = "outcome",                     # для основного уравнения
                                                           # строится предсказание
                     type = "unconditional")               # безусловные предсказания 
wage.star.0 <- wage.star[1]                                # E(y0*)
wage.star.1 <- wage.star[2]                                # E(y1*)

# Рассчитаем оценку условного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*|z)
wage.cond <- predict(model.mle, 
                     newdata = Boris, 
                     part = "outcome",                     # для основного уравнения
                     type = "conditional")                 # условные предсказания 
wage.cond[1]                                               # E(y|z = 0)
wage.cond[2]                                               # E(y|z = 1)
# К сожалению иногда встроенная функция выдает ошибку,
# поэтому приходится осуществлять расчеты вручную
educ.li <- predict(model.mle, 
                   newdata = Boris, 
                   part = "selection",
                   type = "link")
lambda.est.1 <- dnorm(educ.li) / pnorm(educ.li)            # оценка отношения Миллса
lambda.est.0 <- dnorm(educ.li) / pnorm(-educ.li)
wage.cond.1 <- wage.star.1 +
               rho.mle.1 * sigma.mle.1 * lambda.est.1      # E(y*|z = 1)
wage.cond.0 <- wage.star.0 -
               rho.mle.0 * sigma.mle.0 * lambda.est.0      # E(y*|z = 0)

# Посчитаем эффект воздействия
ATET <- wage.cond.1 - wage.cond.0

# Рассчитаем предельный эффект возраста
# на зарплату людей с высшим образованием
# и без него, то есть на E(y*|z=1) и E(y*|z=1)
# соотвественно
coef.s.est <- coef(model.mle)[1:3]
  # аналитически
age.ME <- coef.mle.1["age"] - rho.mle.1 * sigma.mle.1 *
                              (wage.li * lambda.est.1 +
                              lambda.est.1 ^ 2) *
                              coef.s.est ["age"]
