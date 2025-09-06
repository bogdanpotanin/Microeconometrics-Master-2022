# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 1. Модели бинарного выбора
# --------

# Пакеты
library("switchSelection")

# Отключим scientific notation
options(scipen = 999)

# Симулируем данные, содержащую информацию
# о характеристиках заемщика, а также о том,
# наступил ли у него дефолт по ипотеке.
# Все переменные измерены в условных единицах
set.seed(12345)                                            # для воспроизводимости
n     <- 10000                                             # число индивидов в выборке
h     <- data.frame(ind = rep(1:n))                        # датафрейм для хранения данных
h$inc <- runif(n, 0, 1)                                    # доход в условных единицах
h$pay <- runif(n, 0, 1)                                    # ежемесячный платеж
h$age <- runif(n, 0, 1)                                    # возраст
h$ins <- rbinom(n, 1, 0.7)                                 # факт наличия страховки
h$chl <- rbinom(n, 1, 0.6)                                 # факт наличия детей
eps   <- rnorm(n)                                          # случайная ошибка 

beta  <- c(0.6, -3, 2, -5, 3.5, -0.8, 0, 0.8)              # оцениваемые регрессионные 
                                                           # коэффициенты
             

def_li   <- beta[1] +                                      # линейный предиктор (индекс),
            beta[2] * h$inc +                              # отражающий вклад наблюдаемых
            beta[3] * h$pay +                              # факторов в вероятность дефолта
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
# Часть 1. Оценивание пробит модели
#---------------------------------------------------

# Пробит модель
def_formula  <- def ~ inc + pay        +                   # указываем формулу без константы, поскольку
                      age + I(age ^ 2) +                   # она учитывается автоматически
                      ins + chl        +
                      I(chl * inc)
model.probit <- glm(formula = def_formula,   
                    data    = h,                           # датафрейм, из которого берутся 
                                                           # зависимая и независимые переменные
                    family  = binomial(link = "probit"))   # тип оцениваемой бинарной регрессии: в 
                                                           # данном случае пробит
summary(model.probit)                                      # визуализация результата

# Информация об оценках
coef.probit   <- coef(model.probit)                        # оценки регрессионных коэффициентов
as.cov.probit <- vcov(model.probit)                        # оценка асимптотической
                                                           # ковариационной матрицы
                                                           # оценок коэффициентов
as.sd.probit  <- sqrt(diag(as.cov.probit))                 # оценки асимптотических
                                                           # стандартных отклонений
                                                           # оценок коэффициентов

# Сравним истинные значения коэффициента
# и соответствующие оценки
data.frame("True" = beta, "Estimates" = coef.probit)

# Оценим асимптотические доверительные
# интервалы для регрессионных коэффициентов
  # с помощью функции
confint(model.probit,                                      # модель
        level = 0.9)                                       # уровень доверия
  # вручную
cbind("5 %"  = coef.probit - qnorm(0.95) * as.sd.probit,
      "95 %" = coef.probit + qnorm(0.95) * as.sd.probit)

# Протестируем гипотезы о равенстве
# нулю коэффициентов вручную
t       <- coef.probit / as.sd.probit                      # тестовая статистика
p.value <- 2 * pmin(pnorm(t), 1 - pnorm(t))                # p-value теста
p.value <- round(p.value, 5)

# Агрегируем результаты
cbind("Estimate"   = coef.probit,
      "Std. Error" = as.sd.probit,
      "z value"    = t, 
      "Pr(>|z|)"   = p.value)

#---------------------------------------------------
# Часть 2. Оценивание вероятностей в пробит модели
#---------------------------------------------------

# Получим оценки (предсказанные значения) вероятностей
# дефолта для каждого индивида в выборке
probs <- predict(model.probit,                             # модель как объект, возвращаемый glm()
                 type = "response")                        # тип предсказания: в данном
head(probs, 5)                                             # случае вероятность

# Оценим линейный предиктор (индекс)
li <- predict(model.probit,
              type = "link")

# Рассчитаем оценки вероятностей при помощи
# оценок линейных индексов вручную
probs <- pnorm(li)
head(probs, 5)

# Оценивание вероятности дефолта для
# конкретного индивида
Boris <- data.frame(inc = 0.2,                             # характеристики Бориса
                    pay = 0.1,
                    age = 0.3,
                    ins = 1,
                    chl = 1)
prob.Boris <- predict(model.probit,                        # объект возвращаемый glm()
                      newdata = Boris,                     # датафрейм, с использованием значений
                                                           # которого будут предсказаны вероятности
                      type    = "response")                # предсказываем вероятности

# Проверим гипотезу о вероятности дофолта Бориса
# H0: P(У Бориса наступил дефолт) = 0.1
H0_prob       <- 0.1
prob.Boris.se <- predict(model.probit,                     # достанем стандартную ошибку, то есть
                         newdata = Boris,                  # оценку асимптотического стандартного
                         type    = "response",             # отклонения вероятности занятости Бориса
                         se.fit  = TRUE)$se.fit            # применяя аргумент se.fit       
t             <- (prob.Boris - H0_prob) / prob.Boris.se    # считаем статистку теста
p_value       <- 2 * min(pnorm(t), 1 - pnorm(t))           # p-value теста

# Построим симметричный асимптотический 90%-й
# доверительный интервал для вероятности занятости Бориса
ci <- prob.Boris + qnorm(0.95) * prob.Boris.se * c(-1, 1)

# Оценим асимптотическую дисперсию вероятности
# дефолта Бориса вручную, используя дельта метод
Boris$def       <- 1                                       # из технических соображений присваиваем
                                                           # единичное значение зависимой переменной
x.Boris         <- model.frame(formula(model.probit),      # формируем вектор регрессоров для Бориса
                               data = Boris)
x.Boris         <- as.numeric(x.Boris)
li.Boris        <- sum(x.Boris * coef.probit)              # оцениваем линейный индекс Бориса
prob.Boris.grad <- dnorm(li.Boris) * x.Boris               # оцениваем градиент вероятности дефолта
                                                           # Бориса по регрессионным коэффициентам
prob.Boris.se   <- sqrt(prob.Boris.grad %*%                # используем дельта-метод
                        as.cov.probit   %*% 
                        prob.Boris.grad)

# Рассмотрим более общий способ тестирования гипотез

# Оценим пробит модель с помощью альтернативной функции
model.probit2 <- msel(formula  = def_formula,              # формула   
                      data     = h,                        # данные
                      cov_type = "hessian")                # способ оценивания
                                                           # асимптотической
                                                           # ковариационной матрицы

# Запрограммируем функцию, о равенстве 
# нулю которой тестируется гипотеза
# H0: P(У Бориса наступил дефолт) = prob_H0
# H0: P(У Бориса наступил дефолт) - prob_H0 = 0
test_fn <- function(object, prob_H0 = 0.5, newdata = NULL)
{
  prob <- predict(object,                                 # модель
                  newdata = newdata,                      # данные
                  type    = "prob",                       # оцениваем вероятность 
                                                          # того, что
                  group   = 1)                            # P(Y = 1 | X = newdata)
  val        <- prob - prob_H0
  return(val)
}

# Применим тест
test <- test_msel(model.probit2, 
                  fn      = test_fn, 
                  fn_args = list(prob_H0 = 0.1,
                                 newdata = Boris))
summary(test)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 2.1.    Рассчитайте вероятность дефолта для
#         индивида с медианными характеристиками
# 2.2*.   Повторите предыдущий пункт не используя
#         функцию predict()
# 2.3.    Оцените асимптотическую дисперсию оценки вероятности
#         дефолта для индивида с медианными характеристиками,
#         используя дельта метод
#         1*)    используйте аналитическую формулу градиента
#         2**)   используйте численное дифференцирование

#---------------------------------------------------
# Часть 3. Оценивание предельных эффектов
#          в пробит модели
#---------------------------------------------------

library("margins")                                         # библиотека для расчета
                                                           # предельных эффектов
                                                           # Предельные эффекты

# Для корректного оценивания предельных эффектов
# с учетом переменных взаимодействия их следует
# задавать без I(). Для нелинейных эффекто,
# напротив, необходимо указывать I()
model.probit <- glm(formula = def ~ inc + pay        +
                                    age + I(age ^ 2) +     # нелинейный эффект 
                                    ins + chl        +
                                    chl * inc,             # переменная взаимодействия 
                    data    = h,
                    family  = binomial(link = "probit"))
coef.probit  <- coef(model.probit)


probit.me    <- margins(model     = model.probit,          # объект, возвращаемый glm(), в котором   
                                                           # хранятся полученные нами ранее оценки
                        variables = NULL,                  # переменная, по которой считается
                                                           # предельный эффект: можно указать вектор
                                                           # переменных или поставить NULL (по умолчанию), 
                                                           # чтобы получить предельные эффекты сразу по 
                                                           # всем независимым переменным
                        type      = "response")            # на что строится предельный эффект:
                                                           # вероятность или линейный индекс

# Посмотрим оценки предельных эффектов на 
# вероятность дефолта
head(probit.me$dydx_inc, 5)                                # предельный эффект по доходу
head(probit.me$dydx_age, 5)                                # предельный эффект по возрасту
# Визуализируем распределение
# предельных эффектов по дозоду
hist(probit.me$dydx_inc, 
     main = "Распределение предельных эффектов
             дохода на вероятность дефолта",
     xlab = "Предельный эффект",
     ylab = "Доход",
     col  = "cornsilk1")

# По поводу средних предельных эффектов можно
# получить следующие результаты:
# factor - переменная, предельный эффект которой
#          на вероятность подлежит рассмотрению
# AME    - средний предельный эффект
# SE     - оценка асимптотического стандартного
#          отклонения среднего предельного эффекта
# p      - p-value теста о равенстве предельного
#          эффекта нулю
# level  - скольки процентный асимптотический
#          доверительный интервал должен быть 
#          построен
# lower  - нижняя граница (level * 100) процентного
#          асимптотического доверительно интервала
#          для предельного эффекта 
# upper  - верхняя граница (level * 100) процентного
#          асимптотического доверительно интервала
#          для предельного эффекта 
summary(probit.me,                                         # объект, возвращаемый функцией 
                                                           # margins(), в котором
                                                           # хранятся полученные нами 
                                                           # ранее оценки
        level = 0.9)                                       # скольки процентный 
                                                           # асимптотический доверительный 
                                                           # интервал для предельного эффекта 
                                                           # будет построен

# Оценим предельные эффекты для Бориса
me.Boris <- margins(model = model.probit,                  # характеристики Бориса
                    at    = Boris,                         # в форме датафрейма
                    type  = "response")

# По непрерывным переменным
summary(me.Boris)                                          # предельные эффекты 
                                                           # для бинарных переменных
                                                           # в выдаче указаны неправильно
# Рассчитаем предельных эффект по факту
# наличия страховки вручную
# P(Дефолт | X, ins = 1) - P(Дефолт | X, ins = 0)
Boris.noins      <- Boris                                  # создаем индивида идентичного Борису,
Boris.noins$ins  <- 0                                      # но без страховки, а затем
prob.Boris.noins <- predict(model.probit,                  # считаем его вероятность дефолта
                            newdata = Boris.noins, 
                            type    = "response")
me.ins           <- prob.Boris - prob.Boris.noins          # рассчитываем предельный эффект
print(me.ins)

# Автоматизированный расчет предельного эффекта
# факта наличия страховки и тестирование гипотезы
# о равенстве этого предельного эффекта нулю
me.ins2 <- predict(model.probit2,                          # модель
                   newdata = Boris,                        # данные
                   me      = "inc",                        # переменная, по которой 
                                                           # считается предельный 
                                                           # эффект изменения
                   eps     = c(0, 1),                      # от 0 до 1
                   group   = 1,                            # вероятности P(Y = 1 | X)
                   type    = "prob", 
                   test    = TRUE)                         # о предельном эффекте 
                                                           # тестируется гипотеза
                                                           # о равенстве нулю
summary(me.ins2)

# Посчитаем некоторые предельные эффекты вручную
li.Boris <- predict(model.probit,                          # оценка линейного индекса
                    newdata = Boris,                       # для Бориса
                    type    = "link")
den      <- dnorm(li.Boris)                                # функция плотности в точке
                                                           # оценки линейного индекса
me.pay   <- den * coef.probit["pay"]                       # предельный эффект по платежам
me.inc   <- den * (coef.probit["inc"] +                    # предельный эффект по доходу
                   Boris$chl * coef.probit["inc:chl"])
me.age   <- den * (coef.probit["age"] +                    # предельный эффект по возрасту
                   2 * Boris$age * coef.probit["I(age^2)"])

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 3.1.    Для индивида со средними характеристиками
#         оцените предельный эффект возраста на
#         вероятность дефолта
# 3.2*.   Повторите предыдущий пункт не используя
#         функцию margins()
# 3.3**.  Протестируйте значимость рассчитанного вами
#         предельного эффекта не используя margins(),
#         воспользовавшись дельта-методом

#---------------------------------------------------
# Часть 4. Логистическая регрессия и 
#          отношение шансов
#---------------------------------------------------

# Оценим логистическую регрессию
model.logit <- glm(formula = def_formula,                  # для краткости возьмем формулу
                   data    = h,                            # прямиком из оцененной ранее
                                                           # пробит модели
                   family  = binomial(link = "logit"))         
summary(model.logit)

# Достанем оценки коэффициентов
coef.logit <- coef(model.logit)

# Отношение шансов = (вероятность успеха) / (вероятность неудачи)
#                    P(y = 1)             / P(y = 0)                     
step <- 0.1                                                # приращение для переменных

# Изменение отношения шансов
# [P(y = 1 | xk + step) / P(y = 0 | xk + step)] /
# [P(y = 1 | xk)        / P(y = 0 | xk)]

# Оценим, во сколько раз, при прочих равных,
# изменится отношение шансов при
OR.pay <- exp(coef.logit["pay"] * step)                    # изменении платежей на step
OR.ins <- exp(coef.logit["ins"])                           # получении страховки
OR.age <- exp(coef.logit["age"] * step +                   # изменении возраста на step
              coef.logit["I(age^2)"] * step ^ 2 +
              2 * coef.logit["I(age^2)"] * 
              model.logit$data$age * step)
head(OR.age, 5)

# Тестирование гипотез об отношениях шансов

# Оценим модель альтернативным способом
model.logit2 <- msel(formula  = def_formula,
                     data     = h,
                     cov_type = "hessian",
                     marginal = list(logistic = 0))

# Функция для расчета средних отношений шансов
fn_odds <- function(object, step = 1, variable, is_mean = TRUE)
{
  # Вероятности без приращения
  prob <- predict(object, type = "prob", group = 1)
  
  # Вероятности с приращением
  newdata             <- object$data
  newdata[, variable] <- newdata[, variable] + step
  prob_step           <- predict(object, 
                                 type    = "prob", 
                                 group   = 1, 
                                 newdata = newdata)
  
  # Отношения шансов
  odds <- (prob_step / (1 - prob_step)) / (prob / (1 - prob))
  val  <- odds - 1
  
  # Результат
  if (is_mean)
  {
    return(mean(na.omit(val)))
  }
  return(val)
}

# Расчет отношений шансов и тестирование гипотез
odds <- test_msel(model.logit2, 
                  fn      = fn_odds,
                  fn_args =  list(step     = 0.1, 
                                  variable = "pay",
                                  is_mean  = TRUE))
summary(odds)

# Отношения шансов для возраста в отношении каждого наблюдения
odds2 <- test_msel(model.logit2, 
                   fn      = fn_odds,
                   fn_args =  list(step     = 0.1, 
                                   variable = "age",
                                   is_mean  = FALSE))
summary(odds2)

# 4.1.    Оцените, во сколько раз, при прочих равных,
#         изменится отношение шансов вероятности дефолта
#         при изменении
#         1)    переменной pay на 0.2
#         2*)   переменной inc на 0.1
#         3*)   переменной age на 0.1
# 4.2**.  Используя дельта-метод постройте асимптотический
#         доверительный интервалы для изменения отношений шансов
#         вследствие изменения переменной pay на 0.1

#---------------------------------------------------
# Часть 5. Сравнение моделей
#---------------------------------------------------

library("DescTools")                                       # расчет псевдо
                                                           # коэффициентов
                                                           # детерминации
library("pROC")                                            # ROC-curve

# Сравним логит и пробит модели
# разными способами

# 1. Информационные критерии
rbind(probit = AIC(model.probit),                          # сравнение по AIC
      logit  = AIC(model.logit))
rbind(probit = BIC(model.probit),                          # сравнение по BIC
      logit  = BIC(model.logit))

# 2. Псевдо коэффициенты детерминации
rbind(probit = PseudoR2(model.probit, which = "all"),
      logit  = PseudoR2(model.logit,  which = "all"))

# 3. Прогностическая точность
threshold <- 0.5                                           # пороговая вероятность
  # пробит модель
probs.probit <- predict(model.probit,  type = "response")  # вероятность
def.probit   <- as.numeric(probs.probit >= threshold)      # дефолт
tail(cbind(true       = h$def,                             # сравнение предсказаний
           prediction = def.probit), 8)                    # и истинных значений
correct.probit <- mean(def.probit == h$def)                # доля верных прогнозов
  # логит модель
probs.logit   <- predict(model.logit,  type = "response")  # вероятность
def.logit     <- as.numeric(probs.logit >= threshold)      # дефолт
correct.logit <- mean(def.logit == h$def)                  # доля верных прогнозов
  # наивная модель
def.p         <- mean(h$def)                               # доля дефолтов
correct.naive <- max(def.p, 1 - def.p)                     # доля верных прогнозов
  # сравнение (доля верных прогнозов в процентах)
rbind(probit = correct.probit,                             # пробит модель
      logit  = correct.logit,                              # логит модель
      naive  = correct.naive) * 100                        # наивная модель

# Сравнение прогностической точности с 
# использованием тестовой и тренировочной выборок

# Разделим выборку на тестовую и тренировочную
set.seed(777)                                              # для воспроизводимости
n         <- nrow(h)                                       # число наблюдений
n.train   <- round(0.7 * n)                                # число наблюдений в
                                                           # тренировочной и
n.test    <- n - n.train                                   # тестовой выборках
ind.train <- sample(1:n, n.train)                          # индексы наблюдений, попадающих
                                                           # в тренировочную выборку
h.train   <- h[ind.train, ]                                # тренировочная выборка
h.test    <- h[-ind.train, ]                               # тестовая выборка

# Оценим на тренировочных выборках
  # пробит модель
model.probit <- glm(formula = formula(model.probit),      
                    data    = h.train,
                    family  = binomial(link = "probit"))
  # логит модель
model.logit <- glm(formula = formula(model.probit),
                   data    = h.train,
                   family  = binomial(link = "logit"))

# Сравним прогностическую точность
# на тестовых выборках
threshold <- 0.5                                           # порог отсечения
  # пробит модель
probs.probit   <- predict(model.probit, 
                          newdata = h.test,
                          type    = "response")            # вероятность
def.probit     <- as.numeric(probs.probit >= threshold)    # дефолт
correct.probit <- mean(def.probit == h.test$def)           # доля верных прогнозов
  # логит модель
probs.logit    <- predict(model.logit, 
                          newdata = h.test, 
                          type    = "response")            # вероятность
                                                      
def.logit      <- as.numeric(probs.logit >= threshold)     # дефолт
correct.logit  <- mean(def.logit == h.test$def)            # доля верных прогнозов
  # наивная модель
def.p          <- mean(h.test$def)                         # доля дефолтов
correct.naive  <- max(def.p, 1 - def.p)                    # доля верных прогнозов
  # сравнение (доля верных прогнозов в процентах)
rbind(probit = correct.probit,                             # пробит модель
      logit  = correct.logit,                              # логит модель
      naive  = correct.naive) * 100                        # наивная модель

# Чувствительности (true positive rate) и 
# специфичности (true negative rate):
# Чувствительность -   доля верных прогнозов среди тех, у
#                      кого истинное значение бинарной
#                      зависимой переменной равняется 1
# Специфичность    -   доля верных прогнозов среди тех, у
#                      кого истинное значение бинарной
#                      зависимой переменной равняется 0

# Функция для оценки точности предсказаний
acc <- function(x_true, x_pred)
{
  pred_true <- x_true == x_pred                                  # проверяем верность
  # каждого предсказания
  
  df <- data.frame(correct     = mean(pred_true, na.rm = TRUE),  # доля верных прогнозов
                   sensitivity = mean(pred_true[x_true == 1]),   # чувствительность
                   specificity = mean(pred_true[x_true == 0]))   # специфичность
  
  return(df)
}

# Оценивание
acc(h.test$def, def.probit)                                      # пробит модель
acc(h.test$def, def.logit)                                       # логит модель

# ROC-curve
roc(h.test$def ~ probs.probit,                                   # пробит модель
    plot = TRUE, print.auc = TRUE, lwd = 5)
roc(h.test$def ~ probs.logit,                                    # логит модель
    plot = TRUE, print.auc = TRUE, lwd = 5)

#---------------------------------------------------
# Часть 6. Линейная вероятностная модель
#---------------------------------------------------

# Подключим необходиму библиотеку
library("lmtest")

# Оценим параметры линейной вероятностной модели
model.lm <- lm(formula = def_formula,
               data    = h)
summary(model.lm)

# Предскажем вероятности
probs <- predict(model.lm)
plot(sort(probs),                                          # график оценок вероятностей
     xlab = "Наблюдение",                           
     ylab = "Оценка вероятности дефолта",
     col  = "darkcyan")                                  
lines(rep(0, length(probs)), lwd = 3, col = "brown3")

# Построим 90% бутстрапированный доверительный 
# интервал для регрессионных коэффициентов
boot.iter           <- 100                                 # число бутстрап итераций
coef.lm             <- coef(model.lm)                      # МНК оценки
coef.boot           <- matrix(NA,                          # матрица, в которую будем 
                              nrow = boot.iter,            # сохранять оценки 
                              ncol = length(coef.lm))      # коэффициентов
colnames(coef.boot) <- names(coef.lm)
head(coef.boot)

for(i in 1:boot.iter)                       
{
  h.rows.ind      <- sample(1:n, n, replace = TRUE)        # случайным образом отбираем
                                                           # строки датафрейма h
  h.boot          <- h[h.rows.ind, ]                       # формируем новый датафрейм
                                                           # из выбранных строк h
  model.boot      <- lm(formula = formula(model.lm),       # оцениваем модель по      
                        data    = h.boot)                  # новой выборке
                                  
  coef.boot[i, ]  <- coef(model.boot)                      # получаем оценки коэффициентов
}
head(coef.boot)
# Рассмотрим бутстрапированный ДИ 
# для коэффициента при 'inc'
q.L.boot <- quantile(coef.boot[, "inc"], 0.05)             # левая граница
q.R.boot <- quantile(coef.boot[, "inc"], 0.95)             # правая граница
ci.boot  <- c(q.L.boot, q.R.boot)                          # объединяем результат
print(ci.boot)

# Оценим ковариационную матрицу
 # с помощью бутстрапа
cov.boot <- cov(coef.boot)
sd.boot  <- sqrt(diag(cov.boot))
  # аналитически с поправкой на гетероскедастичность
X              <- as.matrix(model.frame(model.lm))
X[, 1]         <- 1
colnames(X)[1] <- "(Intercept)"
p.est          <- X %*% coef.lm
cov.lm.1       <- 0
for (i in 1:nrow(X))
{
  cov.lm.1 <- cov.lm.1 + X[i, ] %*% t(X[i, ]) * p.est[i] * (1 - p.est[i])
}
cov.lm.2         <- qr.solve(t(X) %*% X)
cov.lm           <- cov.lm.2 %*% cov.lm.1 %*% cov.lm.2
colnames(cov.lm) <- rownames(cov.lm)
sd.lm            <- sqrt(diag(cov.lm))
  # сравнение стандартных ошибок
cbind(bootstrap = sd.boot,
      adjust    = sd.lm)
  # протестируем гипотезы о коэффициентах с корректной
  # асимптотической ковариационной матрицей
coeftest(x = model.lm, vcov = cov.boot)

# Автоматизированный расчет асимптотической ковариационной
# матрицей с помощью сэндвича
model.lm2 <- msel(formula2 = def_formula, data = h)
summary(model.lm2)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 6.1.    Посчитайте предельные эффекты по каждой
#         из переменных на вероятность занятости
# 6.2*.   Оцените дисперсии случайных ошибок и
#         рассчитайте долю отрицательных значений
#         среди них
# 6.3.    Постройте 95% бутстрапированный ДИ для коэффициента
#         при переменной pay
# 6.4*.   Постройте 95% бутстрапированный ДИ для
#         вероятности того, что Борис испытает дефолт

#---------------------------------------------------
# Часть 7. Тестирование гипотез о
#          нескольких параметрах
#---------------------------------------------------

library("lmtest")                                          # дополнительные функции для 
                                                           # тестирования гипотез

# Оценим полную модель
model.probit.F <- glm(formula = def ~ inc + pay        +
                                      age + I(age ^ 2) +
                                      ins + chl        +
                                      chl * inc,
                      data    = h,                                  
                      family  = binomial(link = "probit")) 
summary(model.probit.F)

# Пример №1. Проверим гипотезу:
# H0: beta age     = 0
#          age ^ 2 = 0
model.probit.R <- glm(formula = def ~ inc + pay +          # ограниченная модель
                                      ins + chl +          # без переменных
                                      chl * inc,           # на возраст
                      data    = h,                                  
                      family  = binomial(link = "probit")) 
summary(model.probit.R)

# Проверим гипотезу при помощи разных тестов,
# где p-value указан как "Pr(>Chisq)"
lrtest(model.probit.F, model.probit.R)                     # тест отношения 
                                                           # правдоподобия                        
waldtest(model.probit.F, model.probit.R, test = "Chisq")   # теста Вальда

# Пример №2. Проверим гипотезу:
# H0: beta age     = -5
#          age ^ 2 = 3

# Оценим ограниченную модель, учитывая ограничения за 
# счет использования фиксированных значений коэффициентов
# для возраста и возраста в квадрате, с помощью аргумента 
# offset в функции glm()
model.probit.R <- glm(formula = def ~ inc + pay +          # ограниченная модель
                                      ins + chl +
                                      chl * inc,
                      offset  = I(-5 * age) +              # фиксированная часть 
                                I(3 * age ^ 2),            # формулы, где 
                                                           # на коэффициенты 
                                                           # накладывается задаваемые
                                                           # нулевой гипотезой
                                                           # ограничения на параметры
                      data    = h,                                  
                      family  = binomial(link = "probit")) 
summary(model.probit.R)

# Проверим гипотезу LR тестом
lrtest(model.probit.F, model.probit.R)

# Пример №3. Проверим гипотезу: 
# H0: beta inc = -2 * pay
#          age = -5

model.probit.R <- glm(formula = def ~ I(inc - pay / 2) +   # объединяем две 
                                                           # переменные в одну      
                                      ins + chl        +
                                      I(age ^ 2)       +
                                      chl * inc,
                      offset  = I(-5 * age),               # фиксированная часть 
                                                           # формулы, где 
                                                           # на коэффициенты 
                                                           # накладывается задаваемые
                                                           # нулевой гипотезой
                                                           # ограничения на параметры
                      data    = h,                                  
                      family  = binomial(link = "probit")) 
summary(model.probit.R)

# Проверим гипотезу LR тестом, поскольку функция
# waldtest() не учитывает ограничения, введенные 
# за счет аргумента offset в фукнции glm()
lrtest(model.probit.F, model.probit.R)

#---------------------------------------------------
# Часть 8. Тестирование гипотез о возможности
#          оценивания общей модели для групп
#---------------------------------------------------

# Проверим, можно ли оценивать совместную модели
# для людей с детьми и без них
# H0: коэффициенты в моделях для людей с детьми
#     и без - не различаются

# Оценим ограниченную модель
model.probit.R <- glm(formula = def ~ inc + pay        +
                                      age + I(age ^ 2) +
                                      ins + chl,
                      data    = h,
                      family  = binomial(link = "probit"))

# Оценим полную модель как комбинацию двух моделей
  # модель только для людей с детьми
model.probit.F1 <- glm(formula = def ~ inc + pay        +
                                       age + I(age ^ 2) +
                                       ins,
                       data    = h[h$chl == 1, ],
                       family  = binomial(link = "probit"))
  # модель только для людей без детей
model.probit.F0 <- glm(formula = def ~ inc + pay        +
                                       age + I(age ^ 2) +
                                       ins,
                       data    = h[h$chl == 0, ],
                       family  = binomial(link = "probit"))

# Считаем логарифмы правдоподобия 
# полной и ограниченной моделей
lnL.F <- logLik(model.probit.F1) + logLik(model.probit.F0) # логарифм правдоподобия 
                                                           # полной модели
lnL.R <- logLik(model.probit.R)                            # логарифм правдоподобия 
                                                           # ограниченной модели

# Тестируем гипотезу
r       <- length(model.probit.R$coefficients) - 2         # число ограничений
t       <- 2 * (lnL.F - lnL.R)                             # статистика теста
p.value <- as.numeric(1 - pchisq(t, df = r))               # p-value теста

# Альтернатинвый вариант тестирования
model.probit.F <- glm(formula = def ~ inc + pay        +
                                      age + I(age ^ 2) +
                                      ins + chl        +
                                      I(inc * chl)     + 
                                      I(pay * chl)     +
                                      I(age * chl)     +
                                      I(age ^ 2 * chl) +
                                      I(ins * chl),
                      data    = h,
                      family  = binomial(link = "probit"))
  # вручную
lnL.F   <- logLik(model.probit.F)
r       <- length(model.probit.R$coefficients) - 2         # число ограничений
t       <- 2 * (lnL.F - lnL.R)                             # статистика теста
p.value <- as.numeric(1 - pchisq(t, df = r))               # p-value теста
  # с помощью встроенной функции
lrtest(model.probit.F, model.probit.R)

#---------------------------------------------------
# Часть 9. Сравнение коэффициентов между моделями
#---------------------------------------------------

# Важно - для работы методов из данного раздела необходимо скачать версию
#         пакета switchSelection не менее 2.1.0 по ссылке:
#         https://github.com/bogdanpotanin/switchSelection

# Формула
def_formula2 <- def ~ inc + pay        +
                      age + I(age ^ 2) +
                      ins

# Добавляем идентифицирующую переменную в данные
h$crossind <- 1:nrow(h)

# Модель для людей с детьми
model1 <- msel(def_formula2,
               data = h[h$chl == 1, ])

# Модель для людей без детей
  model0 <- msel(def_formula2,
               data = h[h$chl == 0, ])

# Функция для сравнения коэффициентов
fn_chl <- function(object)
{
  coef1 <- coef(object[[1]], eq = 1)
  coef0 <- coef(object[[2]], eq = 1)
  diff  <- coef1 - coef0
  return(diff)
}

# Тестирование гипотезы о равенсвте коэффициентов между моделями
test_chl <- test_msel(list(model1, model0), fn = fn_chl)
summary(test_chl)

# Функция для сравнения средних предельных эффектов
fn_chl2 <- function(object, variable)
{
  me1  <- predict(object[[1]], me = variable, type = "prob", group = 1)
  me0  <- predict(object[[2]], me = variable, type = "prob", group = 1)
  diff <- mean(na.omit(me1)) - mean(na.omit(me0))
  return(diff)
}

# Тестирование гипотезы о равенсвте средних предельных
# эффектов между моделями
test_chl2 <- test_msel(list(model1, model0), 
                       fn      = fn_chl2, 
                       fn_args = list(variable = "pay"))
summary(test_chl2)

#---------------------------------------------------
# Часть 10. Самостоятельное оценивание пробит модели
#---------------------------------------------------

# Будем оценивать пробит регрессию не
# предполагая заранее наличие константы
ProbitLnL <- function(beta, y, X)                          # функция правдоподобия
{
  beta          <- matrix(beta, ncol = 1)                  # вектор оптимизируемых параметров, то есть
                                                           # регрессионных коэффициентов,
                                                           # переводим в матрицу с одним столбцом
  y_est         <- X %*% beta                              # оценка математического ожидания 
                                                           # латентной переменной
  
  n_obs         <- nrow(X)                                 # количество наблюдений
  
  L_vec         <- matrix(NA, nrow = n_obs,                # вектор столбец вкладов наблюдений
                          ncol = 1)                        # в функцию правдоподобия
  
  is_y_0        <- y == 0                                  # вектор условий y = 0
  is_y_1        <- y == 1                                  # вектор условий y = 1
  
  L_vec[is_y_1] <- pnorm(y_est[is_y_1])                    # вклад наблюдений для которых yi = 1
  L_vec[is_y_0] <- 1 - pnorm(y_est[is_y_0])                # вклад наблюдений для которых yi = 0
  
  lnL           <- sum(log(L_vec))                         # логарифм функции правдоподобия
  
  return(lnL)
}

# Пробит регрессия
Probit <- function(formula,                                # формула
                   data)                                   # датафрейм                
{
  d       <- model.frame(formula, data)                    # извлекаем переменные согласно формуле
  y       <- as.matrix(d[, 1], ncol = 1)                   # зависимая переменная как первый
                                                           # столбец в d
  X       <- as.matrix(d[, -1])                            # независимые переменные как все переменные
                                                           # из data кроме зависимой
  X       <- cbind(1, X)                                   # добавим константу
  
  x0_n    <- ncol(X)                                       # число оцениваемых параметров
  
  X_names <- names(d)[-1]                                  # имена независимых переменных
  X_names <- c("Intercept", X_names)                       # добавляем имя константе
  
  result  <- optim(par     = rep(0, x0_n),                 # в качестве начальных точек возьмем нули
                   method  = "BFGS",                       # численный метод оптимизации
                   fn      = ProbitLnL,                    # максимизируемая функция правдоподобия
                   control = list(maxit   = 10000,         # чтобы минимизационную задачу превратить
                                  fnscale = -1,            # в максимизационную умножаем функцию на -1
                                  reltol  = 1e-10),        # установим достаточно высокую точность          
                   hessian = TRUE,                         # вернем Гессиан функции
                   X = X, y = y)                           # аргументы оптимизируемой функции 
  
  
  beta.est        <- result$par                            # оценки коэффициентов
  names(beta.est) <- X_names                               # сопоставляем имена для оценок коэффициентов
  
  as.cov.est           <- qr.solve(-result$hessian)           # оценка асимптотической ковариационной
  colnames(as.cov.est) <- X_names                          # матрицы полученных оценок
  rownames(as.cov.est) <- X_names                          # сопоставляем имена
  
  return_list <- list("beta" = beta.est,                   # возвращаем оценки коэффициентов и
                      "cov" = as.cov.est,                  # асимптотической ковариационной матрицы
                      "data" = data,                       # возвращаем использовавшийся датафрейм
                      "lnL" = result$value,                # возвращаем логарифм функции правдоподобия
                      "X" = X,                             # возвращаем матрицу регрессоров
                      "y" = y,                             # возвращаем зависимую переменную     
                      "formula" = formula)                 # возвращаем использовавшуюся формулу                          
  
  class(return_list) <- "probit"                           # для удобства назначим класс
                                                           # возвращаемой из функции переменной
  return(return_list)                                      # возвращаем результат                               
}

# Воспользуемся созданной функцией
model        <- Probit(def ~ inc + pay + 
                             age + I(age ^ 2) +
                             ins + chl + I(chl * inc),
                data = h)                                    
beta.est     <- model$beta                                 # получаем оценки коэффициентов
beta.cov.est <- model$cov                                  # получаем оценку асимптотической
                                                           # ковариационной матрицы оценок
                                                           # коэффициентов

# Сравним полученные оценки с теми, что возвращает
# встроенная функция
data.frame("glm"    = coef.probit,
           "Probit" = beta.est)

# Метод, красиво обобщающий результаты
# для пробит регрессии
summary.probit <- function(model)                          # пишем перегрузку метода summary()
{                                                          # для объектов класса probit
  n <- length(model$z)                                     # число наблюдений
  
  cat("Probit model estimation results\n")                 # заголовок
  
  cat(paste("Observations:", n, "\n"))                     # пишем число наблюдений
  cat(paste("Log-likelihood:",                             # пишем логарифм функции
            round(model$lnL, 3), "\n"))                    # правдоподобия
  
  cat("---\n")                                             # для красоты
  
  cat("Coefficients:\n")                                   # для красоты
  
  beta.est   <- model$beta                                 # достаем оценки коэффициентов
  
  as_std_est <- sqrt(diag(model$cov))                      # вычисляем оценки асимптотических
                                                           # стандартных отклонений
                                                           # оценок beta
  
  z       <- beta.est / as_std_est                         # считаем тестовые статистики
  p_value <- 2 * pmin(pnorm(z), 1 - pnorm(z))              # рассчитываем p-value
  
  m       <- length(beta.est)                              # число оцениваемых коэффициентов
  
  stars   <- rep("", m)                                    # звездочки для красоты
                                                     
  stars[p_value <= 0.001]                     <- "***"
  stars[(p_value > 0.001) & (p_value < 0.01)] <- "**"
  stars[(p_value > 0.01) & (p_value <= 0.05)] <- "*"
  stars[(p_value > 0.05) & (p_value <= 0.1)]  <- "."
  
  df <- data.frame("Estimate" = beta.est,                  # складываем информацию в датафрейм
                   "As.Std" = as_std_est,
                   "Test.Statistic" = z,
                   "p-value" = p_value,
                   "Significance" = stars)
  is_df_num_col       <- sapply(df, is.numeric)            # отмечаем numeric столбцы датафрейма и
  df[, is_df_num_col] <- round(df[, is_df_num_col], 5)     # округляем их до 5-го знака после точки
  print(df)                                                # печатаем датафрейм
  cat("---\n")                                             # для красоты
  cat(paste("Signif. codes:  0 ‘***’ 0.001 ‘**’",          # для понятности  
            "0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n"))                      
}
# Воспользуемся созданной функцией для
# получения наглядного результата
summary(model)

# Любой предельный эффект нетрудно рассчитать
# с помощью численного дифференцирования, что очень часто
# бывает удобно использовать в очень сложных моделях,
# где аналитический расчет производной является 
# затруднительным:
ProbitME <- function(model,                                # объект возвращаемый glm()
                     variable,                             # по какой переменной считается
                                                           # предельный эффект
                     newdata = NULL,                       # если NULL, то считаем по исходной
                                                           # выборке, а иначе по newdata датафрейму
                     delta   = NULL)                       # если NULL, то приращение для численного      
                                                           # дифференцирования считается автоматически
{                                                              
  if(is.null(newdata))                                     # если не был подана новый датафрейм,
  {                                                        # то считаем предельные вероятности
    newdata <- model$data                                  # по исходному
  }
  
  f0     <- predict(model,                                 # предсказанные вероятности при исходном
                                                           # значении переменной
                    newdata = newdata, 
                    type    = "response")
  
  my.var <- newdata[, variable]                            # достаем переменную, по которой
                                                           # смотрится предельный эффект
  if (is.null(delta))                                      # выбираем приращение, если
  {                                                        # оно заранее не задано пользователем
    delta <- sqrt(.Machine$double.eps) * my.var            
  }
  
  newdata[, variable] <- my.var + delta                    # добавляем приращение к переменной
                                                           # внутри датафрейма
  f1    <- predict(model,                                  # рассчитываем вероятность успеха
                   newdata = newdata,                      # с учетом приращения в отношении
                   type    = "response")                   # переменной
  
  value <- (f1 - f0) / delta                               # считаем предельный эффект через
                                                           # численную производную методом
                                                           # forward difference
  return(value)                                               
}

# Воспользуемся функцией
model.probit <- glm(formula = def ~ inc + pay +
                              age + I(age ^ 2) +     
                              ins + chl +
                              chl * inc,   
                    data    = h,
                    family  = binomial(link = "probit"))
inc.ME <- ProbitME(model    = model.probit,                # получаем оценки предельных эффектов 
                   variable = "inc")                       # дохода на вероятность дефолта
mean(inc.ME)                                               # оценка среднего предельного эффекта
ProbitME(model    = model.probit,                          # оценка предельного эффекта дохода 
         variable = "inc",                                 # Бориса на вероятность дефолта
         newdata  = Boris)
me.Boris$dydx_inc[1]                                       # сравним с полученным ранее с
                                                           # помощью команды margins() результатом
# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 10.1.    Измените функцию Probit() таким образом,
#          чтобы константу не нужно было указывать
#          в формуле и она учитывалась автоматически
#          как в функции glm(). Проследите, чтобы функция
#          summary() при этом работала корректно
# 10.2.    Добавьте в функцию Probit() аргумент cov_method,
#          позволяющий:
#          1)    выбирать, какая оценка асимптотической ковариационной 
#                матрицы оценок будет рассчитана: обратный Гессиан, 
#                произведение Якобианов или сэндвич
#          2*)   добавьте возможность получить оценку данной
#                асимптотической ковариационной матрицы с
#                помощью бутстрапа
# 10.3*.   В функции Probit() сделайте так, чтобы
#          константа не оценивалась, а коэффициент при
#          первом из регрессоров равнялся 1. После чего
#          сделайте sigma и mu оцениваемыми параметрами.
#          Укажите, как соотносятся оценки, полученные
#          в изначальной и репараметризированной модели.
#          Примечание: поскольку при первом из регрессоров
#                      коэффициент зафиксирован на 1, то
#                      априорно предполагается, что он
#                      значим и оказывает положительно
#                      вияние на вероятность успеха
# 10.4**.  Запрограммируйте аналитический градиент функции
#          правдоподобия, максимизируемой в рамках пробит
#          модели и подключите его к функции optim() в
#          функции Probit().
# 10.5**.  Напишите модель бинарного выбора, в которой
#          вместо нормального распределения используется
#          распределение Стьюдента с 5 степенями свободы.
#          Предварительно симулируйте процесс генерации
#          данных с соответствующим распределением
#          случайных ошибок. Обеспечьте функцию summary()
# 10.6***. Пусть t1 и t2 - независимые случайные величины, имеющие 
#          распределение Стьюдента с df1 и df2 степенями свободы
#          соответственно. Также, имеется независимая от них бернулиевская
#          случайная величина V, такая, что P(V = 1) = p. Рассмотрим
#          распределение следующей случайной величины:
#          G = V * (t5 - a) + (1 - V) * (t10 + b)
#          Напишите модель бинарного выбора, в которой используется данное
#          распределение и df1, df2, a, b и p являются оцениваемыми параметрами,
#          константа не оценивается и коэффициент при первом регрессоре
#          зафиксирован на единице. Обеспечьте функцию summary()
# 10.7.    Сделайте так, чтобы функция ProbitME позволяла:
#          1*)   Оценивать предельные эффекты для бинарных
#                переменных, то есть принимающих значения
#                0 или 1 как в случае с переменной ins или chl
#          2**)  Оценива предельные эффекты для дамми
#                переменных
#          3)    Использовать central difference метод
#                численного дифференцирования
# 10.8.    Для каждой из реализованных вами ранее моделей
#          бинарного выбора обеспечьте функцию, которая
#          позволяет, по аналогии с margins():
#          1*)   оценивать предельные эффекты
#          2**)  строить доверительные интервалы для
#                предельных эффектов
#          3**)  тестировать гипотезы о равенстве предельных
#                эффектов определенным значениям