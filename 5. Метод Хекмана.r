# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 5. Метод Хекмана
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("mvtnorm")                       # симуляции из многомерного
                                         # нормального распределения

library("numDeriv")                      # численное дифференцирование

library("stringr")                       # работа со строками
library("tm")

library("switchSelection")               # метод Хекмана

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
h$cost[h$cat == 0] <- Inf                            # заменяем ненаблюдаемые
                                                     # значения на бесконечность
model_mle <- mvoprobit(                              
  formula = cat ~ age + hobby,                       # уравнение отбора
  formula2 = cost ~ hobby + income,                  # основное уравнение
  data = h,                                          # данные
  estimator = "ml")                                  # метод расчета ММП
summary(model_mle)                                   # результат оценивания
coef_mle <- coef(model_mle,                          # сохраним оценки 
                 type = "coef2", eq2 = 1)            # коэффициентов
cov_mle <- coef(model_mle,                           # оценка ковариации между
                type = "cov12", eq = 1, regime = 0)  # случайными ошибками
sigma_mle <- sigma(model_mle)                        # стандартное отклонение
                                                     # случайной ошибки
rho_mle <- cov_mle / sigma_mle                       # оценка корреляции между
                                                     # случайными ошибками
  # метода Хекмана, основанного на
  # двухшаговой процедуре
model_2st <- mvoprobit(                              
  formula = cat ~ age + hobby,                       # уравнение отбора
  formula2 = cost ~ hobby + income,                  # основное уравнение
  data = h,                                          # данные
  estimator = "2step")                                  # метод расчета ММП
summary(model_2st)                                   # результат оценивания
coef_2st <- coef(model_2st,                          # сохраним оценки 
                 type = "coef2", eq2 = 1)            # коэффициентов
cov_2st <- coef(model_2st,                           # оценка ковариации между
                type = "cov12", eq = 1, regime = 0)  # случайными ошибками
sigma_2st <- sigma(model_2st)                        # стандартное отклонение
                                                     # случайной ошибки
rho_2st <- cov_2st / sigma_2st                       # оценка корреляции между
                                                     # случайными ошибками
                                                     # Сравним оценки и истинные значения
  # регрессионные коэффициенты
data.frame("Real" = beta,                            # истинные значения
           "Least Squares" = coef_ls,                # МНК оценки
           "Heckman MLE" = as.vector(coef_mle),      # оценки ММП Хекмана
           "Heckman 2step" = as.vector(coef_2st))    # оценки двухшагового Хекмана
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
cost_star <- predict(model_mle,                  # модель, оценки которой
                                                 # используются для расчетов
                     newdata = Boris,            # данные
                     group2 = 0,                 # указываем, что работает
                                                 # с уравнением 'y' в режиме 0,
                                                 # смысл режима рассмотрим позже
                     type = "val")               # предсказываем математическое
                                                 # ожидание 'y'

# Рассчитаем оценку условного математического 
# ожидания зависимой переменной основного уравнения,
# то есть E(y*|z)
  # E(y*|z = 0)
cost_cond_1 <- predict(model_mle,                  # модель, оценки которой
                                                   # используются для расчетов
                       newdata = Boris,            # данные
                       group2 = 0,                 # указываем, что работает
                                                   # с уравнением 'y' в режиме 0,
                                                   # смысл режима рассмотрим позже
                       type = "val",               # предсказываем математическое
                       group = 1)                  # накладываем условие z = 1 
  # E(y*|z = 1)
cost_cond_0 <- predict(model_mle,                  # модель, оценки которой
                                                   # используются для расчетов
                       newdata = Boris,            # данные
                       group2 = 0,                 # указываем, что работает
                                                   # с уравнением 'y' в режиме 0,
                                                   # смысл режима рассмотрим позже
                       type = "val",               # предсказываем математическое
                       group = 0)                  # накладываем условие z = 0 
                                                   # ожидание 'y'

# Получим предсказания для уравнения отбора:
  # оценим P(z = 1)
cat_prob <- predict(model_mle, 
                    newdata = Boris, 
                    group = 1,                     # предсказываемое значение
                    type = "prob")                 # предсказываем вероятность
  # оценим линейный индекс без константы
cat_li <- predict(model_mle, 
                  newdata = Boris,                 # для уравнения отбора
                  type = "li")                     # предсказываем линейный индекс
  # добавим константу
cat_li <- cat_li - model_mle$cuts$cat

# Оценим E(y*|z) вручную:
lambda_est_1 <- dnorm(cat_li) / pnorm(cat_li)            # оценка отношения Миллса
lambda_est_0 <- dnorm(cat_li) / pnorm(-cat_li)
cost_star + rho_mle * sigma_mle * lambda_est_1           # E(y*|z = 1)
cost_star - rho_mle * sigma_mle * lambda_est_0           # E(y*|z = 0)

# Рассчитаем предельный эффект увлечения котами
# на затраты на содержание кота на E(y*|z=1)
  # автоматически
me_hobby_1 <- predict(model_mle, group2 = 0, 
                    type = "val", group = 1, 
                    me = "hobby", newdata = Boris, 
                    se = TRUE)
summary(me_hobby_1)
  # аналитически
coef_s_est <- coef(model_mle, eq = 1)
me_hobby_2 <- coef_mle[, "hobby"] - rho_mle * sigma_mle *
                                    (cat_li * lambda_est_1 +
              lambda_est_1 ^ 2) * coef_s_est["hobby"]
  # при помощи численного дифференцирования
eps <- 1e-8
Boris_eps <- Boris
Boris_eps$hobby <- Boris$hobby + eps
cost_cond_1_eps <- predict(model_mle, newdata = Boris_eps, 
                           type = "val", group2 = 0, group = 1)
me_hobby_3 <- (cost_cond_1_eps - cost_cond_1) / eps
  # сравнение
c(auto = me_hobby_1[, "val"], analytical = me_hobby_2, numeric = me_hobby_3)

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

# Вернем привычные значения для затрат на кота
h$cost[is.infinite(h$cost)] <- NA

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
      "mvoprobit" = coef_2st)

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

# Сформируем бинарную зависимую переменную
# на большие затраты на кота
h$cost_b <- h$cost > 2                                    # бинарная переменная
                                                          # на превышающие некоторый
                                                          # объем затраты на кота
summary(h$cost_b)

# Заменим ненаблюдаемые значения на -1
h$cost_b[h$cat == 0] <- -1

# Построим модель
model <- mvoprobit(formula = list(cat ~ age + hobby,
                                  cost_b ~ hobby + income),
                   data = h)
summary(model)

# Предскажем вероятность того, что у индивида
# большие затраты на кота, при условии, что
# у него есть кот
probs <- predict(model, group = c(1, 1), given_ind = 1, type = "prob")
head(probs)

# Посчитаем предельный эффект Хобби на вероятность высоких
# затрат на содержание кота
probs2 <- predict(model, group = c(-1, 1), 
                  type = "prob", me = "hobby")
head(probs2)

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
# Часть 5. Пример оценивания уравнения зарплаты
#---------------------------------------------------

# Загрузим данные
data(cps)
help(cps)
data <- cps

# Оценим уравнение зарплаты
data$lwage[data$work == 0] <- Inf
model <- mvoprobit(formula = work ~ age + nchild + health + bachelor + master,
                   formula2 = lwage ~ age + bachelor + master,
                   data = data)
summary(model)

# Ослабим допущение о нормальном распределении
# случайной ошибки уравнения отбора
data$lwage[data$work == 0] <- Inf
model2 <- mvoprobit(formula = work ~ age + nchild + health + bachelor + master,
                    formula2 = lwage ~ age + bachelor + master,
                    marginal = list(hpa = 3), data = data)
summary(model2)

# Проверим гипотезу о нормальности случайной ошибки
# уравннеия отбора
lrtest(model, model2)
c(parametric = AIC(model), semiparametric = AIC(model2))

# Ослабим допущение о совместном нормальном распределении
# и воспользуемся методом Ньюи
model3 <- mvoprobit(formula = work ~ age + nchild + health + bachelor + master,
                    formula2 = lwage ~ age + bachelor + master,
                    marginal = list(hpa = 3),
                    estimator = "2step", degrees = 2,
                    data = data)
summary(model3)

# Проверим гипотезу о том, что отдача от бакалаврского
# и магистерского уровня образования одинаковая
coef(model, type = "coef2")[[1]][, "bachelor"]
coef(model, type = "coef2")[[1]][, "master"]
fn_test <- function(object)
{
  val <- coef(object, type = "coef2")[[1]][, "master"] -
         coef(object, type = "coef2")[[1]][, "bachelor"] 
  return (val)
}
test <- delta_method(model, fn = fn_test)
summary(test)

#---------------------------------------------------
# Часть 6. Метод Ньюи
#---------------------------------------------------

# Предположим, что распределение случайной ошибки отличается от нормального
set.seed(777)
e1 <- rexp(n, 0.7)
e2 <- runif(n, -2, 2)
e3 <- rt(n, df = 5, ncp = 10)
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
h$cost[h$cat == 0] <- Inf
model_mle <- mvoprobit(formula = cat ~ age + hobby,                     
                       formula2 = cost ~ hobby + income,
                       data = h)
coef_mle <- as.vector(coef(model_mle, 
                           eq2 = 1, type = "coef2"))

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

# Автоматическое оценивание методом Ньюи
model.newey.auto <- mvoprobit

# Сравним оценки
data.frame("Real" = beta,                            # истинные значения
           "Heckman MLE" = coef_mle,                 # оценки ММП Хекмана
           "Newey 1" = coef.newey.1[1:3],            # метод Ньюи 
           "Newey 2" = coef.newey.2[1:3],
           "Newey 3" = coef.newey.3[1:3])    
