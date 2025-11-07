# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 6.1. Модели с усеченной зависимой переменной
# --------

# Для удобства отключим 
# экспоненциальную запись чисел
options(scipen = 999)

# Подключим дополнительные библиотеки
library("crch")                          # регрессия с усечением

library("hpa")                           # моменты усеченного 
                                         # нормального распределения

library("switchSelection")               # данные

# Симулируем данные
library("mvtnorm")

# Симулируем данные
set.seed(777)
n <- 10000
h <- data.frame(inc     = exp(rnorm(n, 10, 0.7)),     # доход
                weather = pmin(rpois(n, 0.5) + 1, 8), # погодные условия (1 - плохие, 8 - лучшие)
                deliv   = round(runif(n, 0, 127)))    # цена доставки

# Случайная ошибка
sigma = sqrt(700000)
epsilon <- rnorm(n, mean = 0, sd = sigma)

# Не усеченная зависимая переменная
beta <- c(500, 0.07, -1000, -30)
h$spend <- beta[1]             + 
           beta[2] * h$inc     +
           beta[3] * h$weather +
           beta[4] * h$deliv   + 
           epsilon

# Точки усечения
tr_left <- 2000

# Убираем усечения
h <- h[h$spend >= tr_left, ]
nrow(h)
head(h)

#---------------------------------------------------
# Часть 1. Оценивание параметров модели
#---------------------------------------------------

# Усеченная регрессия
tr       <- 2000                                     # точка усечения
model_tr <- crch(spend ~ inc + weather + deliv,      # формула
                 data      = h,                      # данные                       
                 left      = tr,                     # нижнее (левое) усечение
                 truncated = TRUE,                   # усеченная регрессия 
                 dist      = "gaussian")             # распределение случайной
                                                     # ошибки
summary(model_tr)                                    # посмотрим результат

# Сохраним оценки модели
est_tr   <- coef(model_tr)                           # достанем оценки
coef_tr  <- est_tr[-length(est_tr)]                  # оценки регрессионных
                                                     # коэффициентов
sigma_tr <- exp(est_tr[length(est_tr)])              # достаем оценку 
                                                     # стандартного отклонения

# Воспользуемся обычным МНК для сравнения
model_lm <- lm(spend ~ inc + weather + deliv,
               data = h)
coef_lm <- coef(model_lm)

# Сравним истинные коэффициенты с оценками
cbind(true = beta, 
      trunc = coef_tr, 
      ls = coef_lm)

#---------------------------------------------------
# Часть 2. Расчет вероятностей и предельных эффектов
#---------------------------------------------------

# Рассчитаем оценку безусловного математического
# расходов для конкретного домохозяйства
Boris <- data.frame(inc     = 60000,
                    weather = 2,
                    deliv   = 1)

# Предскажем ожидания зависимой переменной, 
# то есть E(spend*) = E(xb + e) = xb
spend_est <- predict(model_tr, newdata = Boris)

# Предскажем условное математическое
# ожидание расходов  
  # рассчитаем E(e | spend* >= 2000) = 
  #          = E(e | xb + e >= 2000) =
  #          = E(e | e >= 2000 - xb)
    # с помощью функции
epsilon_E <- truncatedNormalMoment(k = 1,                      # момент
                                   x_lower = tr - spend_est,   # нижнее усечение
                                   x_upper = Inf,              # верхние усечение
                                   mean = 0,                   # математическое ожидание
                                   sd = sigma_tr)              # стандртаное отклонение
    # вручную
a         <- (tr - spend_est) / sigma_tr
lambda    <- dnorm(a) / (1 - pnorm(a))
epsilon_E <- sigma_tr * lambda
  # посчитаем E(spend* | spend* >= 2000)
spend_est_cond <- spend_est + epsilon_E

# Рассмотрим предельные эффекты на
# условное математическое ожидание
ME <- coef_tr * (1 - lambda * (lambda - a))

#---------------------------------------------------
# Часть 3. Учет ошибок спецификации
#---------------------------------------------------

# Загрузим данные
data(cps)
data <- cps

# Помыслим гипотетическую ситуацию, что опрос
# проводился только по высокооплачиваемым индивидам
# и для удобства удалим безработных
data                          <- data[data$work == 1, ]
data$lwage2                   <- data$lwage
tr                            <- 3
data$lwage2[data$lwage2 < tr] <- NA

# Распределение до и после усечения
plot(density(na.omit(data$lwage)))   # до усечения
plot(density(na.omit(data$lwage2)))  # после усечения

# Оценим уравнение зарплаты с помощью 
# метода наименьших квадратов
model1 <- lm(lwage2 ~ age + bachelor + master + health, data = data)
summary(model1)

# Воспользуемся моделью усеченной регрессии
model2 <- crch(lwage2 ~ age + bachelor + master + health, 
               data = data, left = tr, truncated = TRUE)
summary(model2)

# Рассмотрим альтернативное распределение случайных ошибок
model3 <- crch(lwage2 ~ age + bachelor + master + health, 
               data = data, left = tr, truncated = TRUE, 
               dist = "student")
summary(model3)

# Учтем гетероскедастичность по здоровью
model4 <- crch(lwage2 ~ age + bachelor + master + health | health, 
               data = data, left = tr, truncated = TRUE)
summary(model4)

# Учтем и гетероскдедасичность и отклонение от нормальности
model5 <- crch(lwage2 ~ age + bachelor + master + health | health, 
               data = data, left = tr, truncated = TRUE, 
               dist = "student")
summary(model5)

# Сравним качество моделей
data.frame(
  cbind(paste0("model", 2:5),
        AIC = round(c(AIC(model2), AIC(model3), 
                      AIC(model4), AIC(model5)), 4),
        BIC = round(c(BIC(model2), BIC(model3), 
                      BIC(model4), BIC(model5)), 4)))

# Сравнение прогнозной точности

# Разделим выборку на тренеровочную и тестовую
train_ind <- sample(1:nrow(data), size = nrow(data) * 0.7)
train     <- data[train_ind, ]
test      <- data[-train_ind, ]

# Оценим модели
model1.train <- lm(formula(model1), 
                   data = train)
model2.train <- crch(formula(model2), data = train,
                     left = tr, truncated = TRUE)
model3.train <- crch(formula(model3), data = train, 
                     left = tr, truncated = TRUE, dist = "student")
model4.train <- crch(formula(model4), data = train, 
                     left = tr, truncated = TRUE)
model5.train <- crch(formula(model5), data = train, 
                     left = tr, truncated = TRUE, dist = "student")

# Для удобства реализуем функцию для прогнозирования
# условных математических ожиданий
predict.cond <- function(model, newdata)
{
  y.uncond  <- predict(model, newdata = newdata)
  sigma_tr  <- predict(model, type = "scale", newdata = newdata)
  epsilon_E <- truncatedNormalMoment(
    k = 1,                      
    x_lower = (model$cens$left - y.uncond) / sigma_tr,   
    x_upper = (model$cens$right - y.uncond) / sigma_tr) 
  y.cond    <- y.uncond + epsilon_E
  return(y.cond)
}

# Получим прогнозы на тестовой выборке
y1 <- predict(model1.train, newdata = test)
y2 <- predict.cond(model2.train, newdata = test)

# Посчитаем ошибки прогноза по условным
# математическим ожиданиям
e1 <- test$lwage2 - y1
e2 <- test$lwage2 - y2

# Вычислим RMSE
data.frame(rmse1 = sqrt(mean(na.omit(e1) ^ 2)),
           rmse2 = sqrt(mean(na.omit(e2) ^ 2)))

# Сравним точность прогнозов на низкую зарплату,
# попавшую под усечение

# Безулосвные прогнозы
y1 <- predict(model1.train, newdata = test)
y2 <- predict(model2.train, newdata = test)
y3 <- predict(model3.train, newdata = test)
y4 <- predict(model4.train, newdata = test)
y5 <- predict(model5.train, newdata = test)

# Ошибки
e1 <- test$lwage[test$lwage < tr] - y1[test$lwage < tr]
e2 <- test$lwage[test$lwage < tr] - y2[test$lwage < tr]
e3 <- test$lwage[test$lwage < tr] - y3[test$lwage < tr]
e4 <- test$lwage[test$lwage < tr] - y4[test$lwage < tr]
e5 <- test$lwage[test$lwage < tr] - y5[test$lwage < tr]

# RMSE
data.frame(rmse1 = sqrt(mean(na.omit(e1) ^ 2)),
           rmse2 = sqrt(mean(na.omit(e2) ^ 2)),
           rmse3 = sqrt(mean(na.omit(e3) ^ 2)),
           rmse4 = sqrt(mean(na.omit(e4) ^ 2)),
           rmse5 = sqrt(mean(na.omit(e5) ^ 2)))