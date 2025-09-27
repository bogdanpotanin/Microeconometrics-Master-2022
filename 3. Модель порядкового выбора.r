# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 3. Модель порядкового выбора
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("mnorm")                                           # симуляции из многомерного
                                                           # нормального распределения
library("numDeriv")                                        # численное дифференцирование
library("brant")                                           # тест Бранта о parallel lines
library("MASS")                                            # порядковый пробит и логит
library("switchSelection")                                 # порядковый пробит и логит
                                                        
# Воспроизведем процесс генерации данных,
# предполагаемый пробит моделью порядкового выбора

# Симулируем данные
set.seed(123)                                              # для воспроизводимости
n      <- 10000                                            # число наблюдений                  
X      <- rmnorm(n,                                        # симулируем n наблюдений из многомерного
                                                           # нормального распределения
             c(0, 0, 0),                                   # с нулевым вектором математических ожиданий и
             matrix(c(1,    0.2,  0.3,                     # следующей ковариационной матрице
                      0.2,    1, -0.1,
                      0.3, -0.1,    1),
                    ncol = 3,
                    byrow = FALSE))
X[, 2] <- X[, 2] > 0.5                                     # сделаем регрессор X2 бинарным
u      <- rnorm(n, 0, 1)                                   # ошибки из стандартного нормально распределения

# Соберем регрессоры в датафрейм
h <- data.frame("income"     = X[, 1],                     # показатель, отражающий уровень
                                                           # развития профессиональных равыков
                "male"       = X[, 2],                     # пол: мужчина = 1, женщина = 0
                "experience" = X[, 3],                     # показатель, отражающий опыт работы
                "educ"       = runif(n))                   # показатель образования, симулированный из
                                                           # стандартного равномерного распределения
head(h)

# Создадим линейный предиктор, который определяет
# зависимость латентной переменной от различных
# независимых переменных
gamma  <- c(0.1, 0.2, 0.3, -0.05, 0.5, 0)                  # линейный индекс с
z_li   <- gamma[1] * h$income +                            # образованием
          gamma[2] * h$male +                              # полом
          gamma[3] * h$experience +                        # опытом
          gamma[4] * h$experience ^ 2 +                    # квадратом опыта
          gamma[5] * h$income * h$male +                   # взаимодействием навыков и образования
          gamma[6] * h$educ                                # доходом
z_star <- z_li + u

# Введем пороги
mu <- c(0, 0.5, 1.2)

# Создадим наблюдаемую зависимую переменную,
# отражающую оценку индивидом качества
# своего здоровья
z                                       <- rep(0, n)       # наблюдаемое значение переменной
z[(z_star > mu[1]) & (z_star <= mu[2])] <- 1
z[(z_star > mu[2]) & (z_star <= mu[3])] <- 2
z[(z_star > mu[3])]                     <- 3
z        <- matrix(z, ncol = 1)                            # как матрица с одним столбцом 
h$health <- z                                              # добавим в датафрейм переменную 
                                                           # на трудоустройство
summary(as.factor(h$health))                               # посмотрим доли

# Посмотрим на данные
head(h, 5)

#---------------------------------------------------
# Часть 1. Оценивание параметров и тест Бранта
#---------------------------------------------------

# Перекодируем переменную на здоровье в бинарную
h$health_binary <- as.numeric(h$health >= 2)

# Воспользуемся обычно пробит моделью
model_probit <- msel(formula = health_binary ~ income + male +
                                               experience + I(experience ^ 2) +
                                               I(income * male) + educ,
                    data    = h)    
summary(model_probit) 

# Применим порядковую пробит модель
model_oprobit <- msel(formula = health ~ income + male +
                                         experience + I(experience ^ 2) +
                                         I(income * male) + educ, 
                      data    = h)
summary(model_oprobit)

# Сравним истинные и полученные оценки границ
data.frame("True Thresholds" = mu,
           "Estimate"        = model_oprobit$cuts[[1]])

# Сравним точность оценок
gamma_probit  <- coef(model_probit, type = "coef", eq = 1)  # оценки пробит модели без константы
gamma_oprobit <- coef(model_oprobit, type = "coef", eq = 1) # оценки порядковой пробит модели
data.frame("Real"        = gamma,
           "Probit"      = gamma_probit,
           "Ordered"     = gamma_oprobit,
           "MSE Probit"  = (gamma_probit - gamma) ^ 2,
           "MSE Ordered" = (gamma_oprobit - gamma) ^ 2)

# Проведем тест Бранта, предварительно оценив порядковую
# пробит модель с помощью альтернативной функции
model_oprobit2 <- polr(formula = as.factor(health) ~              # зависимую переменную необходимо
                                 income + male +                  # преобразовать в факторную, что
                                 experience + I(experience ^ 2) + # можно сделать прямо в формуле
                                 I(income * male) + educ,                                                 
                      data     = h,                                       
                      method   = "probit")                        # можно заменить на logistic             
brant(model_oprobit2)

# Проведем тест Бранта вручную
h$crossind <- 1:nrow(h)

# Оценим модели
  # первая модель
h$h1 <- as.numeric(h$health >= 1)
m1   <- msel(formula = h1 ~ income + male +
                            experience + I(experience ^ 2) +
                            I(income * male) + educ, 
             data    = h)
  # вторая модель
h$h2 <- as.numeric(h$health >= 2)
m2   <- msel(formula = h2 ~ income + male +
                            experience + I(experience ^ 2) +
                            I(income * male) + educ, 
             data    = h)
# третья модель
h$h3 <- as.numeric(h$health >= 3)
m3   <- msel(formula = h3 ~ income + male +
                            experience + I(experience ^ 2) +
                            I(income * male) + educ, 
             data    = h)

# Функция для тестирования гипотезы о равенстве коэффициентов между моделями
fn_brant <- function(object)
{
  # Коэффициенты моделей
  coef1 <- coef(object[[1]], type = "coef", eq = 1)
  coef2 <- coef(object[[2]], type = "coef", eq = 1)
  coef3 <- coef(object[[3]], type = "coef", eq = 1)
  
  # Разницы коэффициентов
  val <- c(coef1 - coef2, coef2 - coef3)
  
  # Возвращаем разницы
  return(val)
}

# Тестирование гипотезы
test_brant <- test_msel(list(m1, m2, m3), fn = fn_brant, test = "wald")
summary(test_brant)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 1.1.    Разбейте уровни здоровья на две категории
#         альтернативным способом и убедись, что
#         оценки регрессионных коэффициентов окажутся
#         достаточно близки к истине
# 1.2.    Разбейте уровни здоровья на три категории
#         и проверьте, насколько точными окажутся оценки.
#         Сравните результаты данного и предыдущего
#         заданий.
# 1.3*.   Проверьте, как изменятся результаты теста
#         Бранта при:
#         1)     некорректной спецификации
#                линейного индекса
#         2)     гетероскедастичных случайных ошибках
#         3)     если случайная ошибка имеет
#                распределение Стьюдента с пятью
#                степенями свободы
#         4)     если произвести оценивание параметров 
#                модели объединив последние две категории
#                в одну

#---------------------------------------------------
# Часть 2. Расчет вероятностей и предельных эффектов
#---------------------------------------------------

# Оценим вероятности P(health = 2)   
probs_oprobit <- predict(model_oprobit, type = "prob", group = 2)
head(probs_oprobit)

# Посчитаем вероятность попадания не менее, чем во вторую группу P(health >= 2)
probs <- cbind(Probit  = as.numeric(predict(model_probit, type = "prob", 
                                            group = 1)),
               Ordered = as.numeric(predict(model_oprobit, type = "prob", # P(health >= 2) =
                                            group = 2) +                  # P(health = 2)  +
                                    predict(model_oprobit, type = "prob", # P(heath  = 3)
                                            group = 3)))
head(probs)

# Рассчитаем условную вероятность P(health = 2 для конкретного индивида
Boris   <- data.frame(income = 1, male = 1, experience = 1, educ = 0.5)
p_Boris <- predict(model_oprobit, type = "prob", group = 2, newdata = Boris)

# Оценим предельный эффект опыта на P(health = 2)  
me_experience <- predict(model_oprobit, type = "prob", group = 2, 
                         me = "experience")
me_Boris      <- predict(model_oprobit, type = "prob", group = 2, 
                         me = "experience", newdata = Boris)

# Протестируем значимость среднего предельного эффекта опыта на P(health = 2)  
ame_experience <- predict(model_oprobit, type = "prob", group = 2, 
                          me = "experience", test = mean)
summary(ame_experience)

# Оценим средний предельный эффект пола на P(health = 2)  
ame_male <- predict(model_oprobit, type = "prob", group = 2, 
                    me = "male", test = mean, eps = c(0, 1))
summary(ame_male)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 2.1.    Оцените предельный эффект:
#         1)     пола
#         2)     дохода
#         3)     образования
# 2.2*.   Повторите предыдущий пункт, используя
#         численную оптимизацию
# 2.3.    Сравните вероятности попадания в последнюю
#         категорию по пробит и порядковой моделям
# 2.4.    Самостоятельно запрограммируйте процедуры,
#         позволяющие оценивать:
#         1**)   порядковую пробит модель с
#                гетероскедастичной случайной ошибкой
#         2***)  систему из двух порядковых уравнений
#         3***)  порядковую модель со случайными эффектами

#---------------------------------------------------
# Часть 3. Порядковая логистическая модель
#---------------------------------------------------

# Применим порядковую логистическую модель
model_ologit <- msel(formula = health ~ income + male +
                                        experience + I(experience ^ 2) +
                                        I(income * male) + educ,
                     data = h, marginal = list(logistic = NULL))
coef_ologit <- coef(model_ologit, type = "coef", eq = 1)
summary(model_ologit)

# При интерпретации отношений шансов учитывайте, что для 
# любого номера категории k изменение в отношениях шансов 
# p(z>k)/p(z<=k) остается прежним. Рассчитаем изменение в 
# отношениях шансов p(z>k)/p(z<=k) при изменении независимой 
# переменной, входящей в основное уравнение, на единицу
OR_educ <- exp(coef_ologit["educ"])

# Тестирование гипотезы о том, что изменение отношений шансов
# для образования значимо отличается от 1
OR_educ_fn <- function(object)
{
  val <- exp(coef(object, eq = 1)["educ"]) - 1
  return(val)
}
test <- test_msel(model_ologit, OR_educ_fn)
summary(test)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 3.1.    Оцените изменение в отношениях шансов при
#         изменении на единицу:
#         1*)    пола
#         2*)    дохода
#         3*)    опыта

#---------------------------------------------------
# Часть 4. Пример на реальных данных
#---------------------------------------------------

# Оценим влияние различных индивидуальных характеристик
# на субъективную оценку здоровья индивида

# Загрузим данные, содержащие информацию о социально-демографических
# характеристиках замужних женщин в возрасте 25-54 лет в 2022 году
data("cps")
help(cps)
data <- cps

# Перекодируем переменную на уровень здоровья таким образом, чтобы
# ее минимальное значение равнялось нулю
data$health <- data$health - min(data$health, na.rm = TRUE)
table(data$health)

# Оценим порядковую пробит модель
model1 <- msel(health ~ age + work + bachelor + master + nchild, data = data)
summary(model1)

# Оценим порядковую логистическую модель
model2 <- msel(health ~ age + work + bachelor + master + nchild,
               data = data, marginal = list(logistic = NULL))
summary(model2)

# Оценим полупараметрическую модель
model3 <- msel(health ~ age + work + bachelor + master + nchild,
               data = data, marginal = list(hpa = 3))
summary(model3)

# Сравним качество моделей
c(probit = AIC(model1), logit = AIC(model2), semiparametric = AIC(model3))

# Посчитаем вероятность среднего уровня здоровья P(health = 2)
prob.1         <- predict(model1, group = 2, type = "prob")
prob.2         <- predict(model2, group = 2, type = "prob")
prob.3         <- predict(model3, group = 2, type = "prob")
prob           <- cbind(prob.1, prob.2, prob.3)
colnames(prob) <- c("probit", "logit", "semiparametric")
head(prob)
plot(x = prob.1, y = prob.3, xlab = "probit", ylab = "semiparametric")

# Посчитаем предельный эффект возраста на вероятность
# среднего уровня здоровья
me.age.1 <- predict(model1, group = 2, type = "prob", me = "age")
me.age.2 <- predict(model2, group = 2, type = "prob", me = "age")
me.age.3 <- predict(model3, group = 2, type = "prob", me = "age")
hist(me.age.1, breaks = 30)
hist(me.age.2, breaks = 30)
hist(me.age.3, breaks = 30)

# Посчитаем средний предельный эффект возраста на вероятность
# среднего уровня здоровья
ame.age.1 <- mean(me.age.1)
ame.age.2 <- mean(me.age.2)
ame.age.3 <- mean(me.age.3)
c(probit = ame.age.1, logit = ame.age.2, semiparametric = ame.age.3)

# Протестируем гипотезу о равенстве нулю соответствующего среднего
# предельного эффекта
me.age.fn <- function(object)
{
  val <- mean(predict(object, group = 2, type = "prob", me = "age"))
  return(val)
}
me.age.test.1 <- test_msel(model1, fn = me.age.fn)
summary(me.age.test.1)
me.age.test.2 <- test_msel(model2, fn = me.age.fn)
summary(me.age.test.2)
me.age.test.3 <- test_msel(model3, fn = me.age.fn)
summary(me.age.test.3)

# Посчитаем средний предельный эффект наличия магистерского образования
# в сравнение с бакалаврским на вероятность отличного здоровья
ame.master.fn <- function(object)
{
  p_master   <- predict(object, group = 4, type = "prob",
                        exogenous = list(bachelor = 0, master = 1))
  p_bachelor <- predict(object, group = 4, type = "prob",
                        exogenous = list(bachelor = 1, master = 0))
  val        <- mean(p_master - p_bachelor)
  return(val)
}
ame.master <- test_msel(model1, ame.master.fn)
summary(ame.master)
