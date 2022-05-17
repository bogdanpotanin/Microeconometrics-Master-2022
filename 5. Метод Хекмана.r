# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 5. Метод Ньюи
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

library("GJRM")                          # бинарный Хекман с копулами

#---------------------------------------------------
# Часть 1. Оценивание параметров
#---------------------------------------------------

set.seed(123)
# Зададим количество наблюдений
n <- 10000
# Истинные значения регрессионных коэффициентов
beta <- c(1, 3, 5)
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
# Симулируем случайные ошибик из двумерного
# нормального распределения
sigma <- 3                                        # стандартное отклонение случайной
                                                  # ошибки основного уравнения
rho <- 0.8                                        # корреляция между случайными ошибками
epsilon_cov <- matrix(c(1, rho * sigma,           # ковариационная матрица 
                        rho * sigma, sigma ^ 2),  # случайных ошибок
                      ncol = 2)
epsilon <- rmvnorm(n, c(0, 0), epsilon_cov)       # случайные ошибки
# Латентная зависимая переменная
  # уравнения отбора 
z_star <- gamma[1] + gamma[2] * X[, 1] + 
          gamma[3] * X[, 2] + epsilon[, 1]
  # основного уравнения
y_star <- beta[1] + beta[2] * X[, 2] + 
          beta[3] * X[, 3] + epsilon[, 2]

# Наблюдаемая зависимая переменная
  # основного уравнения
z <- ifelse(z_star >= 0, 1, 0)
  # уравнения отбора 
y <- ifelse(z == 1, y_star, NA)
# Посмотрим на результат
head(cbind(z, y))

# Сюжет:
# Представим, что вы изучаете, как различные 
# факторы влияют на затраты индивидов на содержание 
# кота. При этом вы наблюдаете объем затрат лишь для
# индивидов, у которых есть кот.
h <- data.frame("cost" = y,                          # затраты на кота
                "cat" = z,                           # факт наличия кота
                "age" = X[, 1],                      # возраст
                "hobby" = X[, 2],                    # увлечение котами
                "income" = X[, 3])                   # доход
head(h)

# Оценим параметры модели при помощи
  # МНК
model_ls <- lm(cost ~ hobby + income,  
               data = h)
summary(model_ls)                                    # результат оценивания
coef_ls <- coef(model_ls)                            # сохраним оценки коэффициентов
sigma(model_ls)                                      # оценка стандартного отклонения
  # метода Хекмана, основанный на ММП
model_mle <- selection(                              
  selection = cat ~ age + hobby,                     # уравнение отбора
  outcome = cost ~ hobby + income,                   # основное уравнение
  data = h,                                          # данные
  method = "ml")                                     # метод расчета ММП
summary(model_mle)                                   # результат оценивания
coef_mle <- coef(model_mle, part = "outcome")        # сохраним оценки коэффициентов
rho_mle <- model_mle$estimate["rho"]                 # оценка корреляции между
                                                     # случайными ошибками
sigma_mle <- model_mle$estimate["sigma"]             # стандартное отклонение
                                                     # случайной ошибки
  # метода Хекмана, основанного на
  # двухшаговой процедуре
model_2st <- selection(                              
  selection = cat ~ age + hobby,                     
  outcome = cost ~ hobby + income,                   
  data = h,                                          
  method = "2step")                                  # метод расчета двухшаговая процедура
summary(model_2st)                                   # результат оценивания
coef_2st <- coef(model_2st, part = "outcome")        # сохраним оценки коэффициентов
coef_2st <- coef_2st[-length(coef_2st)]              # удалим лишний коэффициент
rho_2st <- model_2st$rho                             # оценка корреляции между
                                                     # случайными ошибками
# Сравним оценки и истинные значения
  # регрессионные коэффициенты
data.frame("Real" = beta,                            # истинные значения
           "Least Squares" = coef_ls,                # МНК оценки
           "Heckman MLE" = coef_mle,                 # оценки ММП Хекмана
           "Heckman 2step" = coef_2st)               # оценки двухшагового Хекмана
  # корреляция случайных ошибок
data.frame("Real" = rho,                             # истиннoе значение
           "Least Squares" = 0,                      # МНК оценки
           "Heckman MLE" = rho_mle,                  # оценка ММП Хекмана
           "Heckman 2step" = rho_2st)                # оценка двухшагового Хекмана

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 1.1.    Посмотрите, как изменится результат, если:
#         1)     при корреляции между случайными
#                ошибками, равной 0.1, 0.9, -0.5
#         2)     если не соблюдены exclusion
#                restrictions: в обоих уравнениях
#                одни и те же регрессоры
# 1.2*.   Придумайте и реализуйте собственный
#         пример на симулированных данных
# 1.3.    Проверьте необходимость в оценивании
#         модели с учетом неслучайного отбора
#         при помощи:
#         1) LR теста
#         2*) Wald теста
#         3*) LM теста

#---------------------------------------------------
# Часть 2. Расчет предсказаний и предельных эффектов
#---------------------------------------------------

# Создадим индивида
Boris <- data.frame("cost" = 1,
                    "cat" = 1,
                    "age" = 0.1,
                    "hobby" = 0.3,
                    "income" = 0.5)

# Рассчитаем оценку безусловного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*)
cost_star <- predict(model_mle, 
                     newdata = Boris, 
                     part = "outcome",                   # для основного уравнения
                                                         # строится предсказание
                     type = "unconditional")             # безусловные предсказания   

# Рассчитаем оценку условного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*|z)
cost_cond <- predict(model_mle, 
                     newdata = Boris, 
                     part = "outcome",                   # для основного уравнения
                     type = "conditional")               # условные предсказания 
cost_cond[1]                                             # E(y*|z = 0)
cost_cond[2]                                             # E(y*|z = 1)

# Получим предсказания для уравнения отбора:
  # оценим P(z = 1)
cat_prob <- predict(model_mle, 
                    newdata = Boris, 
                    part = "selection",                  # для уравнения отбора
                    type = "response")                   # предсказываем вероятность
  # оценим линейный индекс
cat_li <- predict(model_mle, 
                  newdata = Boris, 
                  part = "selection",                    # для уравнения отбора
                  type = "link")                         # предсказываем линейный индекс

# Оценим E(y*|z) вручную:
lambda_est_1 <- dnorm(cat_li) / pnorm(cat_li)            # оценка отношения Миллса
lambda_est_2 <- dnorm(cat_li) / pnorm(-cat_li)
cost_star + rho_mle * sigma_mle * lambda_est_1           # E(y*|z = 1)
cost_star - rho_mle * sigma_mle * lambda_est_2           # E(y*|z = 0)

# Рассчитаем предельный эффект увлечения котами
# на затраты на содержание кота на E(y*|z=1)
coef_s_est <- coef(model_mle)[1:3]
  # аналитически
hobby_ME <- coef_mle["hobby"] - rho_mle * sigma_mle *
                                (cat_li * lambda_est_1 +
                                 lambda_est_1 ^ 2) *
                                coef_s_est["hobby"]
  # при помощи численного дифференцирования
eps <- 1e-8
Boris_eps <- Boris
Boris_eps$hobby <- Boris$hobby + eps
cost_cond_eps <- predict(model_mle, 
                         newdata = Boris_eps, 
                         part = "outcome",
                         type = "conditional")
hobby_ME_num <- (cost_cond_eps[2] - 
                 cost_cond[2]) / eps

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 2.1.    Для индивида с произвольными характеристиками
#         оцените предельный эффект возраста на:
#         1)     E(y*)
#         2)     P(z = 1)
#         3)     E(y|z = 1)
#         4)     E(y|z = 0)
#         5*)    P(z = 1|y = 0.5)
# 2.2.    Повторите предыдущее задание добавив
#         в модель:
#         1*)    доход в квадрате
#         2*)    взаимодействие между возрастом
#                и доходом

#---------------------------------------------------
# Часть 3. Дхвухшаговая процедура
#---------------------------------------------------

# Реализуем самостоятельно двухшаговую процедуру

# На первом шаге
model_1 <- glm(cat ~ age + hobby,                       # оценим пробит модель
               data = h,
               family = binomial(link = "probit"))
cat_li <- predict(model_1, type = "link")               # оценим линейный индекс
h$lambda <- dnorm(cat_li) / pnorm(cat_li)               # оценим лямбду

  # На втором шаге
model_2 <- lm(cost ~ hobby + income + lambda,           # используем МНК включив
              data = h)                                 # оценку лямбды как регрессор                           
coef_2 <- coef(model_2)                                 # достаем оценки коэффициентов

# Сравним полученные результаты со
# встроенной функцией
rbind("Our" = coef_2[1:length(coef_2st)], 
      "sampleSelection" = coef_2st)

# Оценим стандартное отклонение случайной ошибки
res <- model_2$residuals                                # достаем вектор остатков
delta <- h$lambda * (h$lambda + cat_li)                 # вспомогательное выражение
sigma_est <- sqrt(mean(res ^ 2) +                       # оцениваем стандартное отклонение
                  mean(delta[h$cat == 1]) *             # случайной ошибки
                  coef_2["lambda"] ^ 2)
names(sigma_est) <- "sigma"
print(sigma_est)

# Оценим корреляцию между случайными ошибками
rho_est <- coef_2["lambda"] / sigma_est

# Оценим ковариационную матрицу с учетом
# гетероскедастичности
Delta <- diag(delta[h$cat == 1])                        # вспомогательные выражения
Delta_adj <- diag(1 - (rho ^ 2) * delta[h$cat == 1])
X_star <- cbind(1, h$hobby[h$cat == 1],                 # регрессоры основного
                   h$income[h$cat == 1],                # уравнения, включая
                   h$lambda[h$cat == 1])                # лямбду и константу
W <- cbind(1, h$age[h$cat == 1],                        # регрессоры уравнения отбора,
              h$hobby[h$cat == 1])                      # включая константу
gamma_cov <- vcov(model_1)                              # ковариационная матрица оценок
                                                        # коэффициентов уравнения отбора
M <- (t(X_star) %*% Delta %*% W)                        # вспомогательные выражения
Q <- (rho_est ^ 2) * (M %*% gamma_cov %*% t(M))
X_prod_inv <- solve(t(X_star) %*% X_star)
beta_cov <- sigma_est ^ 2 *                             # оценка ковариационной матрицы
            X_prod_inv %*%                              # оценок коэффициентов основного
            ((t(X_star) %*% Delta_adj                   # уравнения
              %*% X_star) + Q) %*% 
            X_prod_inv
beta_std <- sqrt(diag(beta_cov))                        # оценки стандартных отклонений
# Запишем регрессоры и их стандартные отклонения
cbind(coef_2, beta_std)

# Получим красивую выдачу с учетом коррекции
# ковариационной матрицы
library("lmtest")
coeftest(model_2, vcov = beta_cov)

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 3.1.    Проверьте значимость регрессоров
# 3.2.    Проверьте гипотезу о возможности оценивания
#         МНК модели вместо двухшаговой процедуры

#---------------------------------------------------
# Часть 4. Метод Хекмана с бинарной
#          зависимой переменной
#---------------------------------------------------
h$cost_b <- h$cost > 2                                    # бинарная переменная
                                                          # на превышающие некоторый
                                                          # объем затраты на кота
summary(h$cost_b)

# Построим модель
model <- gjrm(formula = list(cat ~ age + hobby,
                             cost_b ~ hobby + income),
              Model = "BSS",
              BivD = "N",
              margins = c("probit", "probit"),
              data = h)
summary(model)

# Проверим значимость корреляции
rho_est <- model$theta                                  # оценка корреляции между случайными ошибками
data.frame("Rho real" = rho,                            # сравним истинное значение корреляции
           "Rho estimate" = rho_est)                    # с её оценкой
cov_est <- solve(model$fit$hessian)                     # оценка асимптотической ковариационной матрицы
std_rho <- sqrt(cov_est["theta.star", "theta.star"])    # оценка стандартной ошибки оценки корреляции
p_value_rho <- 2 * min(pnorm(rho_est / std_rho),        # p-value теста о равенстве корреляции между
                       1- pnorm(rho_est / std_rho))     # случайными ошибками нулю
cov_u_est <- matrix(c(1, rho_est,                       # оценка ковариационной матрицы
                      rho_est, 1),                      # совместного распределения случайных ошибок
                    ncol = 2)

# Таблицу с копулами и их описанием можно найти в:
# Marra G, Wyszynski K (2016), Semi-Parametric Copula Sample Selection Models for 
# Count Responses. Computational Statistics and Data Analysis, 104, 110-129.
# Другие работы авторов, связанные с пакетом GJRM
citation("GJRM")

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# Выполните задания для индивида с произвольными характеристиками: 
# 4.1.    Используя функцию predict() оцените вероятность
#         того, что затраты индивида на содержание кота
#         превысят пороговое значение
# 4.2*.   Оцените вероятность того, что у индивида есть
#         кот и затраты на его содержание превысят
#         пороговую величину
# 4.3*.   Оцените вероятность того, что затраты на содержание
#         кота превысят пороговую величину, при условии,
#         что у него есть кот

#---------------------------------------------------
# Часть 5. Метод Ньюи
#---------------------------------------------------

# Напишем вспомогательную функции для подсчета MSE 
# при leave-one-out кросс валидации
loocv <- function(fit)
{
  h <- lm.influence(fit)$h
  loocv_value <- mean((residuals(fit)/(1 - h)) ^ 2)
  
  return(loocv_value)
}

# Предположим, что распределение случайной ошибки отличается от нормального
set.seed(777)
e1 <- rt(n, df = 5)
e2 <- rt(n, df = 5, ncp = 10)
e3 <- rt(n, df = 5)
b1 <- rbinom(n, 1, 0.5)
epsilon <- cbind((b1 * (e1 - 5) + (1 - b1) * (e2 + 5)) / 11.3, 
                 (b1 * (e1 - 5) + (1 - b1) * (e3 + 5)) / 1.7)
plot(density(epsilon[, 1]))
plot(density(epsilon[, 2]))
cor(epsilon)

# Латентная зависимая переменная
  # уравнения отбора 
z_star <- gamma[1] + gamma[2] * X[, 1] + 
          gamma[3] * X[, 2] + epsilon[, 1]
  # основного уравнения
y_star <- beta[1] + beta[2] * X[, 2] + 
          beta[3] * X[, 3] + epsilon[, 2]

# Наблюдаемая зависимая переменная
  # основного уравнения
h$cat <- ifelse(z_star >= 0, 1, 0)
  # уравнения отбора 
h$cost <- ifelse(h$cat == 1, y_star, NA)

# Попробуем воспользоваться обычным методом Хекмана
model_mle <- selection(                              
  selection = cat ~ age + hobby,                     # уравнение отбора
  outcome = cost ~ hobby + income,                   # основное уравнение
  data = h,                                          # данные
  method = "ml")                                     # метод расчета ММП
summary(model_mle)                                   # результат оценивания
coef_mle <- coef(model_mle, part = "outcome")        # сохраним оценки коэффициентов
rho_mle <- model_mle$estimate["rho"]                 # оценка корреляции между
                                                     # случайными ошибками
sigma_mle <- model_mle$estimate["sigma"]             # стандартное отклонение
                                                     # случайной ошибки

# Сравним истинные значения с оценками
data.frame("Real" = beta,                            # истинные значения
           "Heckman MLE" = coef_mle)                 # оценки ММП Хекмана

# Применим метод Ньию, для быстроты используя
# на первом шаге обычную пробит модель
model.probit <- glm(formula = cat ~ age + hobby,                                                 
                   data = h,                                         
                   family = binomial(link = "probit")) 

# Оценим линейный индекс и сгладим его
# отношением Миллса
li <- predict(model.probit)
h$lambda <- dnorm(li) / pnorm(li)

# Применим метод Ньюи с различными степенями
  # 1
model.newey.1 <- lm(cost ~ hobby + income + 
                           lambda, data = h)
coef.newey.1 <- coef(model.newey.1)
mse.1 <- loocv(model.newey.1)
  # 2
model.newey.2 <- lm(cost ~ hobby + income + 
                           lambda + I(lambda ^ 2), data = h)
coef.newey.2 <- coef(model.newey.2)
mse.2 <- loocv(model.newey.2)
  # 3
model.newey.3 <- lm(cost ~ hobby + income + 
                           lambda + I(lambda ^ 2) + I(lambda ^ 3), data = h)
coef.newey.3 <- coef(model.newey.3)
mse.3 <- loocv(model.newey.3)

# Найдем модель с наименьшим MSE
which.min(c(mse.1, mse.2, mse.3))

# Сравним оценки
data.frame("Real" = beta,                            # истинные значения
           "Heckman MLE" = coef_mle,                 # оценки ММП Хекмана
           "Newey 1" = coef.newey.1[1:3],            # метод Ньюи 
           "Newey 2" = coef.newey.2[1:3],
           "Newey 3" = coef.newey.3[1:3])                            