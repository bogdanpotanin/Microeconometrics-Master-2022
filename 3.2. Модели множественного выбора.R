# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 3.2. Модели множественного выбора
# --------

# Для удобства отключим 
# экспоненциальную запись чисел
options(scipen = 999)

library("mlogit")                                          # мультиномиальный логит
library("switchSelection")                                 # мультиномиальная пробит
library("stringr")                                         # работа со строками
library("tm")                                              # обработка данных
library("mvtnorm")
library("EnvStats")


# Симулируем данные
set.seed(123)                                              # для воспроизводимости
n   <- 10000                                               # число наблюдений  
h   <- data.frame(income = exp(rnorm(n, 10, 0.7)),         # показатель дохода
                  health = pmin(rpois(n, 3), 5),           # показатель здоровья
                  age    = round(runif(n, 20, 100)))      
eps <- cbind(revd(n), revd(n), revd(n))                    # случайные ошибки                                    

# Соберем в датафрейм регрессоры, различные
# для всех альтернатив
h$price.Car      <- exp(rnorm(n, 1, 0.5))
h$comfort.Car    <- pmin(rpois(n, 3), 5)
h$price.Taxi     <- exp(rnorm(n, 2, 0.3))
h$comfort.Taxi   <- pmin(rpois(n, 3.5), 5)
h$price.Public   <- exp(rnorm(n, 0.5, 0.1))
h$comfort.Public <- pmin(rpois(n, 1), 5)
head(h, 5)

# Создадим несколько латентных переменных, каждая
# из которых отражает предпочтения в отношении того
# или иного вида транспорта
beta_Car     <- c(0.1, 0.000025, 0.3, 0.01)
beta_Taxi    <- c(0.2, 0.000015, 0.2, 0.015)
beta_Public  <- c(3, -0.00002, 0.5, -0.02)
beta_d       <- c(-0.3, 0.3)
y_li_Car     <- beta_Car[1] +                              # линейный предиктор Машины
                h$income * beta_Car[2] +
                h$health * beta_Car[3] +
                h$age * beta_Car[4] +
                h$price.Car * beta_d[1] +
                h$comfort.Car * beta_d[2]
y_star_Car   <- y_li_Car + eps[, 1]                        # латентная переменная Машины
y_li_Taxi    <- beta_Taxi[1] +                             # линейный предиктор Такси
                h$income * beta_Taxi[2] +
                h$health * beta_Taxi[3] +
                h$age * beta_Taxi[4] +
                h$price.Taxi * beta_d[1] +
                h$comfort.Taxi * beta_d[2]
y_star_Taxi  <- y_li_Taxi + eps[, 2]                       # латентная переменная Такси
y_li_Public  <- beta_Public[1] +                           # линейный предиктор
                h$income * beta_Public[2] +                # общественного транспорта
                h$health * beta_Public[3] +
                h$age * beta_Public[4] +
                h$price.Public * beta_d[1] +
                h$comfort.Public * beta_d[2]
y_star_Public <- y_li_Public + eps[, 3]                    # латентная переменная 
                                                           # общественного транспорта

# Сформируем зависимую переменную
h$transport[(y_star_Car    >= y_star_Taxi) &               # те, кто выбрал Машину
            (y_star_Car    >= y_star_Public)] <- "Car"
h$transport[(y_star_Taxi   >= y_star_Car) & 
            (y_star_Taxi   >= y_star_Public)] <- "Taxi"    # те, кто выбрал Такси
h$transport[(y_star_Public >= y_star_Car) & 
            (y_star_Public >= y_star_Taxi)]   <- "Public"  # те, кто выбрал
                                                           # общественный транспорт
summary(as.factor(h$transport))

# Посмотрим на данные
head(h, 7)

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 1.1*.   Повторите приведенный пример,
#         рассмотрев собственное множество
#         из пяти альтернатив
# 1.2**.  Симулируйте случайные ошибки из совместного
#         нормального распределения и оцените
#         устойчивость результата

#---------------------------------------------------
# Часть 1. Оценивание параметров модели
#---------------------------------------------------

# Подготовим данные
h1 <- dfidx(h,                                   # исходный датафрейм
            varying = 4:9,                       # индексы регрессоров, которые
                                                 # разнятся между альтернативами
                                                 # имена этих регрессоров должны 
                                                 # иметь формат имя.альтернатива
            shape   = "wide",                 
            choice  = "transport")               # переменная, отражающая 
                                                 # выбранную альтернативу

# Оценим параметры модели
model_cmlogit <- mlogit(transport ~              # формула, включающая
                        price + comfort |        # различающиеся и
                        income + health + age,   # общие регрессоры
                        data = h1)               # преобразованный датафрейм
summary(model_cmlogit)   

#---------------------------------------------------
# Часть 2. Оценивание вероятностей и
#          предельных эффектов
#---------------------------------------------------

# Оценим вероятности по всей выборке
probs_cmlogit <- predict(model_cmlogit,          # без аргумента newdata результат
                         newdata = h1)           # может оказаться неправильным
head(probs_cmlogit, 5)

# Оценим вероятности для индивида
Boris <- data.frame(income         = 10000,      # характеристика Бориса
                    health         = 5,
                    age            = 25,
                    price.Car      = 3,
                    comfort.Car    = 3,
                    price.Taxi     = 6,
                    comfort.Taxi   = 5,
                    price.Public   = 1,
                    comfort.Public = 1,
                    transport      = "Car")      # указываем любое из возможных
                                                 # значений зависимой переменной,
                                                 # что не скажется на результате
Boris1 <- dfidx(Boris, varying = 4:9,            # перекодируем Бориса в
                shape = "wide", "transport")     # нужный формат
probs_Boris <- predict(model_cmlogit,            # предскажем вероятности
                       newdata = Boris1)

# Оценим предельный эффект на вероятность
# выбора Бориса для различных переменным
  # по доходу
ME_income <- effects(model_cmlogit,
                     covariate = "income",
                     data      = Boris1)
  # по цене
ME_price <- effects(model_cmlogit,
                    covariate = "price",
                    data      = Boris1)

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 2.1.    Для индивида с произвольными характеристиками
#         оцените:
#         1)     вероятность того, что он поедет на Такси
#         2)     отношение вероятностей поездки на 
#                Такси и на Машине
# 2.2.    Проверьте гипотезу о возможности исключить
#         из модели:
#         1)     Переменную на здоровье
#         2)     Все специфические для альтернатив регрессоры
#         3)     Все регрессоры, не различающиеся
#                между альтернативами
#         Подсказка: воспользуйтесь LR тестом
# 2.3.    Для индивида с произвольными характеристиками
#         оцените предельный эффект здоровья на:
#         1)     вероятность поездки на Такси
#         2*)    отношение вероятностей поездки на
#                Такси и на Машине
# 2.4.    Повторите предыдущее задание рассмотрев
#         предельный эффект комфорта
# 2.5.    Для индивида с произвольными характеристиками 
#         проверьте гипотезу о том, что он с равной
#         вероятностью:
#         1**)   Поедет на Такси или на Машине
#         2**)   Поедет на любом из видов транспорта
# 2.6.    Постройте 95% доверительный интервал для
#         предельного эффекта:
#         1**)  цены поездки на Такси на вероятность
#               поездки на Машине
#         2**)  цены поездки на Такси на вероятность
#               поездки на Такси
#         3**)  здоровья индивид на вероятность
#               поездки на Такси

#---------------------------------------------------
# Часть 3. Тестирование гипотезы о независимости
#          от посторонних альтернатив IIA
#---------------------------------------------------

# Тестирование гипотезы об IIA
mlogit::hmftest(model_cmlogit,                   # независимость выбора между
                z = c("Taxi", "Car"))            # Такси и Машиной транспортом
                                                 # от возможности использовать 
                                                 # Общественным транспортом

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 3.1*.   Повторите пример с пятью альтернативами и
#         проверьте возможность удалить любые две 
#         из них

#---------------------------------------------------
# Часть 4. Множественная пробит модель
#---------------------------------------------------

# Очень долго оцениваем параметры модели
model_cmprobit <- mlogit(transport ~ price + comfort |       # синтаксис такой же, как    
                                     income + health + age,  # для логит модели
                         data = h1, 
                         probit = TRUE)                      # добавляем этот аргумент            
summary(model_cmprobit)

# В случае, если отсутствуют специфические для альтернатив
# характеристики, оценивание можно осуществить альтернативным
# образом
h$transport1                        <- 0
h$transport1[h$transport == "Taxi"] <- 1
h$transport1[h$transport == "Car"]  <- 2
model_mnprobit <- msel(formula3 = transport1 ~ income + health + age, 
                       type3 = "probit", data = h)
summary(model_mnprobit)

# ЗАДАНИЯ (* - средне, ** - сложно, *** - очень сложно)
# 4.1***. Дождитесь окончания расчетов пробит модели
#         множественного выбора

#---------------------------------------------------
# Часть 5. Дополнительный пример со сравнением
#          моделей порядкового и множественного
#          выбора
#---------------------------------------------------

# Подключим данные
data("cps")
help(cps)
data <- cps

# Распределение индивидов по уровням образования
table(data$educ)

# Зададим формулу
frm <- educ ~ age + I(age ^ 2 / 100) + sbachelor + smaster

# Построим порядковую логит модель выбора уровня образования
model.ologit <- msel(formula = frm, data = data, marginal = list(logistic = 0))
summary(model.ologit)

# Построим порядковую пробит модель выбора уровня образования
model.oprobit <- msel(formula = frm, data = data)
summary(model.oprobit)

# Построим мультиномиальную логит модель выбора уровня образования
model.mnlogit <- msel(formula3 = frm, data = data, type3 = "logit")
summary(model.mnlogit)

# Построим мультиномиальную пробит модель выбора уровня образования
model.mnprobit <- msel(formula3 = frm, data = data, type3 = "probit")
summary(model.mnprobit)

# Сравним модели по информационным критериям
AIC(model.ologit, model.oprobit, model.mnlogit, model.mnprobit)
BIC(model.ologit, model.oprobit, model.mnlogit, model.mnprobit)

# Оценим вероятность наличия магистерского образования
prob.ologit.2   <- predict(model.ologit,   group  = 2, type = "prob")
prob.oprobit.2  <- predict(model.oprobit,  group  = 2, type = "prob")
prob.mnlogit.2  <- predict(model.mnlogit,  group3 = 2, type = "prob_mn")
prob.mnprobit.2 <- predict(model.mnprobit, group3 = 2, type = "prob_mn")
head(cbind(prob.ologit.2, prob.oprobit.2, prob.mnlogit.2, prob.mnprobit.2))

# Оценим предельные эффекты возраста на вероятность магистерского образования
me.age.ologit   <- predict(model.ologit,     group = 2, 
                           type = "prob",    me    = "age")
me.age.oprobit  <- predict(model.oprobit,    group = 2, 
                           type = "prob",    me    = "age")
me.age.mnlogit  <- predict(model.mnlogit,    group3 = 2, 
                           type = "prob_mn", me     = "age")
me.age.mnprobit <- predict(model.mnprobit,   group3 = 2, 
                           type = "prob_mn", me     = "age")
head(cbind(me.age.ologit, me.age.oprobit, me.age.mnlogit, me.age.mnprobit))

# Протестируем значимость среднего предельного эффекта
# возраста на вероятность бакалаврского образования
  # модели порядкового выбора
ame.fn.ordered <- function(object)
{
  me.age <- predict(object, group = 2, type = "prob", me = "age")
  return(mean(me.age))
}
test.ologit  <- test_msel(model.ologit,  fn = ame.fn.ordered)
test.oprobit <- test_msel(model.oprobit, fn = ame.fn.ordered)
  # модели множественного выбора
ame.fn.multimomial <- function(object)
{
  me.age <- predict(object, group3 = 2, type = "prob_mn", me = "age")
  return(mean(me.age))
}
test.mnlogit  <- test_msel(model.mnlogit,  fn = ame.fn.multimomial)
test.mnprobit <- test_msel(model.mnprobit, fn = ame.fn.multimomial)
  # результаты
summary(test.ologit)
summary(test.oprobit)
summary(test.mnlogit)
summary(test.mnprobit)

# Повторим тестирование с помощью бутстрапа
# для порядковой пробит модели
  # для воспроизводимости
set.seed(123)
  # бутстрап
boot.oprobit <- bootstrap_msel(model.oprobit, iter = 100)
head(boot.oprobit$par)
  # тестирование
test2.oprobit <- test_msel(model.oprobit, 
                           fn        = ame.fn.ordered, 
                           test      = "t",
                           bootstrap = boot.oprobit,        
                           method    = "bootstrap",        # способ тестирования
                           se_type   = "bootstrap",        # стандартные ошибки
                           ci        = "percentile",       # доверительный интервал
                           cl        = 0.95)               # уровень доверия
summary(test2.oprobit)
  # воспользуемся альтернативным методом бутстрапа
test3.oprobit <- test_msel(model.oprobit, 
                           fn     = ame.fn.ordered, 
                           test   = "Wald",
                           iter   = 1000,      
                           method = "score")
summary(test3.oprobit)
  # дополнительная информация
help(test_msel)

# Оценим предельные эффекты наличия магистерского уровня образования у супруга
# в сравнении с бакалаврским на вероятность магистерского образования у женщины
me.educ.ologit.fn <- function(object)
{
  p_master   <- predict(object, group = 2, type = "prob", 
                        exogenous = list(sbachelor = 0, smaster = 1))
  p_bachelor <- predict(object, group = 2, type = "prob", 
                      exogenous = list(sbachelor = 1, smaster = 0))
  ame        <- mean(p_master - p_bachelor)
  return (ame)
}
test_educ <- test_msel(model.ologit, fn = me.educ.ologit.fn)
summary(test_educ)

#---------------------------------------------------
# Часть 6. Сравнение численного и аналитического
#          подходов к расчету предельных эффектов
#---------------------------------------------------

# Предельные эффекты удобно рассчитывать численно
# Прсчииаем предельный эффект дохода
# на вероятность каждой из альтернатив
delta               <- sqrt(.Machine$double.eps)        # приращение
Boris1_delta        <- Boris1                           # датафрейм с приращением
Boris1_delta$income <- Boris1_delta$income + delta      # добавляем приращение по
                                                        # доходу в датафрейм
probs_Boris_delta   <- predict(model_cmlogit,           # считаем вероятность с
                               newdata = Boris1_delta)  # учетом приращения
ME_income_num       <- (probs_Boris_delta -             # оцениваем предельный эффект
                        probs_Boris) / delta            # с помощью численного
                                                        # дифференцирования
data.frame("Analytical" = ME_income,                    # сравниваем численный и
           "Numeric"    = ME_income_num)                # аналитический результаты
# Посчитаем предельный эффект цены поездки
# на Машине на вероятность каждой
# из альтернатив
delta         <- sqrt(.Machine$double.eps)              # приращение
Boris1_delta  <- Boris1                                 # датафрейм с приращением
Boris1_is_car <- Boris1_delta$idx$id2 == "Car"          # строки, в которых хранится
# информация о регрессорах,
# специфических для Машины
Boris1_delta$price[Boris1_is_car] <-                    # делаем приращение цены
  Boris1_delta$price[Boris1_is_car] + delta             # поездки на Машине
probs_Boris_delta <- predict(model_cmlogit,             # считаем вероятность с
                             newdata = Boris1_delta)    # учетом приращения
ME_price_num      <- (probs_Boris_delta -               # оцениваем предельный эффект
                      probs_Boris) / delta              # с помощью численного
# дифференцирования
data.frame("Analytical" = ME_price["Car", ],            # сравниваем численный и
           "Numeric"    = ME_price_num)                 # аналитический результаты

# Предельные эффекты можно посчитать и по формуле
# Рассчитаем предельный эффект на вероятность
# использования Такси
# по доходу
coef_cmlogit   <- coef(model_cmlogit)
ME_income_Taxi <- probs_Boris["Taxi"] * 
  (coef_cmlogit["income:Taxi"] - 
     probs_Boris["Car"]    * 0 -
     probs_Boris["Taxi"]   * coef_cmlogit["income:Taxi"] -
     probs_Boris["Public"] * coef_cmlogit["income:Public"])
data.frame("Analytical" = ME_income["Taxi"],                 # сравним результаты
           "Numeric" = ME_income_Taxi)                       # со встроенной функцией
# по цене поездки на 
# общественнрм Транспорте
ME_price_Taxi_Public <- -probs_Boris["Taxi"] *
                         probs_Boris["Public"] *                  
                         coef_cmlogit["price"]
# по цене поездки на Такси
ME_price_Taxi_Taxi <- probs_Boris["Taxi"] * 
                      (1 - probs_Boris["Taxi"]) *
                      coef_cmlogit["price"]
data.frame("Analytical" = ME_price["Taxi", "Taxi"],          # сравним результаты
           "Numeric"    = ME_price_Taxi_Taxi)                # со встроенной функцией