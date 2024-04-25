# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 8. Модель порядкового выбора
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("mvtnorm")                                       # симуляции из многомерного
                                                         # нормального распределения
library("numDeriv")                                      # численное дифференцирование
library("brant")                                         # тест Бранта о parallel lines
library("MASS")                                          # порядковый пробит и логит
library("switchSelection")                               # порядковый пробит и логит
                                                        
# Воспроизведем процесс генерации данных,
# предполагаемый пробит моделью порядкового выбора
# Симулируем данные
set.seed(123)                                            # для воспроизводимости
n <- 10000                                               # число наблюдений                  
X <- rmvnorm(n,                                          # симулируем n наблюдений из многомерного
                                                         # нормального распределения
             c(0, 0, 0),                                 # с нулевым вектором математических ожиданий и
             matrix(c(1, 0.2, 0.3,                       # следующей ковариационной матрице
                      0.2, 1, -0.1,
                      0.3, -0.1, 1),
                    ncol = 3,
                    byrow = FALSE))
X[, 2] <- X[, 2] > 0.5                                   # сделаем регрессор X2 бинарным
X <- cbind(1, X)                                         # добавим константу как дополнительный регрессор
u <- rnorm(n, 0, 1)                                      # ошибки из стандартного нормально распределения

# Соберем регрессоры в датафрейм
h <- data.frame("income" = X[, 2],                       # показатель, отражающий уровень
                                                         # развития профессиональных равыков
                "male" = X[, 3],                         # пол: мужчина = 1, женщина = 0
                "experience" = X[, 4],                   # показатель, отражающий опыт работы
                "educ" = runif(n))                       # показатель веса, симулированный из
                                                         # стандартного равномерного распределения
head(h)

# Создадим линейный индекс, который определяет
# зависимость латентной переменной от различных
# независимых переменных
gamma <- c(0.1, 0.2, 0.3, -0.05, 0.5, 0)                 # линейный индекс с
z_li <- gamma[1] * h$income +                            # образованием
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
z <- rep(0, n)                                           # наблюдаемое значение переменной
z[(z_star > mu[1]) & (z_star <= mu[2])] <- 1
z[(z_star > mu[2]) & (z_star <= mu[3])] <- 2
z[(z_star > mu[3])] <- 3
z <- matrix(z, ncol = 1)                                 # как матрица с одним столбцом 
h$health <- z                                            # добавим в датафрейм переменную 
                                                         # на трудоустройство
summary(as.factor(h$health))                             # посмотрим доли

# Посмотрим на данные
head(h, 5)

#---------------------------------------------------
# Часть 1. Оценивание параметров и тест Бранта
#---------------------------------------------------

# Воспользуемся пробит модель, предварительно перекодировав health
h$health_binary <- as.numeric(h$health >= 2)

model_probit <- glm(formula = health_binary ~ income + male +       # указываем формулу без константы, поскольку
                              experience + I(experience ^ 2) +      # она учитывается автоматически
                              I(income * male) +
                              educ,                                                 
                    data = h,                                       # датафрейм, из которого берутся зависимая
                                                                    # и независимые переменные
                    family = binomial(link = "probit"))             # тип оцениваемой бинарной регрессии: в данном
                                                                    # случае используется пробит регрессия и для ее
                                                                    # замены на логит измените "probit" на "logit"
summary(model_probit) 

# Применим порядковую модель
model_oprobit <- polr(formula = as.factor(health) ~                 # зависимую переменную необходимо
                                income + male +                     # преобразовать в факторную, что
                                experience + I(experience ^ 2) +    # можно сделать прямо в формуле
                                I(income * male) +
                                educ,                                                 
                      data = h,                                       
                      method = "probit")                            # можно заменить на logistic             
summary(model_oprobit)    

# Альтернативный способ c p-value
model_oprobit2 <- mvoprobit(formula = health ~ income + male +
                                      experience + I(experience ^ 2) +
                                      I(income * male) + educ,
                            data = h)
summary(model_oprobit2)

# Сравним истинные и полученные оценки границ
data.frame("True Thresholds" = mu,
           "polr" = model_oprobit$zeta,
           "mvoprobit" = model_oprobit2$cuts[[1]])

# Сравним точность оценок
gamma_probit <- coef(model_probit)[-1]                              # пробит оценки без константы
gamma_oprobit <- coef(model_oprobit)                                # оценки порядкового пробита
data.frame("Real" = gamma,
           "Probit" = gamma_probit,
           "Ordered" = gamma_oprobit,
           "MSE Probit" = (gamma_probit - gamma) ^ 2,
           "MSE Ordered" = (gamma_oprobit - gamma) ^ 2)

# Проведем тест Бранта
brant(model_oprobit)

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

# Оценим вероятности
probs_oprobit <- predict(model_oprobit, type = "probs")             # порядковый пробит
probs_oprobit2 <- predict(model_oprobit2, type = "prob",            # P(health = 2)        
                          group = 2)
probs_probit <- predict(model_probit, type = "response")            # пробит
head(probs_oprobit)

# Сравним вероятности попадания не менее, чем
# во вторую категорию
probs <- data.frame("Probit" = probs_probit,
                    "Ordered" = probs_oprobit[, 3] + probs_oprobit[, 4])
head(probs)

# Оценим предельный эффект опыта на вероятность
# попадания во вторую категорию
  # автоматический способ
ME_experience_1 <- predict(model_oprobit2, type = "prob", 
                           group = 2, me = "experience")
  # ручной расчет
X_oprobit <- model.frame(model_oprobit)[, -1]                       # достаем матрицу независимых переменных
X_oprobit <- as.matrix(X_oprobit)                                   # из порядковой пробит модели
z_li_est <- X_oprobit %*% gamma_oprobit
ME_experience_2 <- (gamma_oprobit["experience"] +                   # часть, обусловленная формой, в которой
                  2 * gamma_oprobit["I(experience^2)"] *            # переменная входит в линейный индекс
                  X_oprobit[, "experience"]) *
                  (dnorm(model_oprobit$zeta[2] - z_li_est) -        # часть, общая для предельных эффектов
                   dnorm(model_oprobit$zeta[3] - z_li_est))         # всех переменных
  # сопоставление
head(cbind(ME_experience_1, ME_experience_2))

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

# Строго говоря в процессе генерации данных нужно заменить
# u <- rnorm(n, 0, 1)
# на
# u <- rlogis(n)

# Применим порядковую модель
model_ologit <- polr(formula = as.factor(h$health) ~                # зависимую переменную необходимо
                               income + male +                      # преобразовать в факторную, что
                               experience + I(experience ^ 2) +     # можно сделать прямо в формуле
                               I(income * male) +
                               educ,                                                 
                      data = h,                                       
                      method = "logistic")                          # можно заменить на logistic             
summary(model_ologit)

# Альтернативный способ c p-value и стандартизацией дисперси
# случайной ошибки к единице
model_ologit2 <- mvoprobit(formula = health ~ income + male +
                                     experience + I(experience ^ 2) +
                                     I(income * male) + educ,
                            data = h, marginal = list(logistic = NULL))
summary(model_ologit2)

# При интерпретации отношений шансов учитывайте, что для 
# любого номера категории k изменение в отношениях шансов 
# p(z>k)/p(z<=k) остается прежним. Рассчитаем изменение в 
# отношениях шансов p(z>k)/p(z<=k) при изменении независимой 
# переменной, входящей в основное уравнение, на единицу
OR_educ <- exp(gamma_oprobit["educ"])

# Тестирование гипотезы о том, что изменение отношений шансов
# для образования значимо отличается от 1
OR_educ_fn <- function(object)
{
  val <- exp(coef(model_ologit2, eq = "health")["educ"]) - 1
  return(val)
}
test <- delta_method(model_ologit2, OR_educ_fn)
summary(test)

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 3.1.    Оцените изменение в отношениях шансов при
#         изменении на единицу:
#         1*)    пола
#         2*)    дохода
#         3*)    опыта

#---------------------------------------------------
# Часть 4. Продвинутый пример
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
model1 <- mvoprobit(health ~ age + work + bachelor + master + nchild, 
                    data = data)
summary(model1)

# Оценим порядковую логистическую модель
model2 <- mvoprobit(health ~ age + work + bachelor + master + nchild,
                    data = data, marginal = list(logistic = NULL))
summary(model2)

# Оценим полупараметрическую модель
model3 <- mvoprobit(health ~ age + work + bachelor + master + nchild,
                    data = data, marginal = list(hpa = 3))
summary(model3)

# Сравним качество моделей
c(probit = AIC(model1), logit = AIC(model2), semiparametric = AIC(model3))

# Посчитаем вероятность среднего уровня здоровья
prob.1 <- predict(model1, group = 2, type = "prob")
prob.2 <- predict(model2, group = 2, type = "prob")
prob.3 <- predict(model3, group = 2, type = "prob")
prob <- cbind(prob.1, prob.2, prob.3)
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
me.age.test.1 <- delta_method(model1, fn = me.age.fn)
summary(me.age.test.1)
me.age.test.2 <- delta_method(model2, fn = me.age.fn)
summary(me.age.test.2)
me.age.test.3 <- delta_method(model3, fn = me.age.fn)
summary(me.age.test.3)

# Посчитаем средний предельный эффект наличия магистерского образования
# в сравнение с базовым на вероятность отличного здоровья
data.master <- data
data.master$bachelor <- 0
data.master$master <- 0
ame.master.1 <- mean(predict(model1, group = 4, type = "prob", 
                             me = "master", eps = c(0, 1),
                             newdata = data.master))
ame.master.2 <- mean(predict(model2, group = 4, type = "prob", 
                             me = "master", eps = c(0, 1),
                             newdata = data.master))
ame.master.3 <- mean(predict(model3, group = 4, type = "prob", 
                             me = "master", eps = c(0, 1),
                             newdata = data.master))
c(probit = ame.master.1, logit = ame.master.2, semiparametric = ame.master.3)

# Посчитаем средний предельный эффект наличия бакалаврского образования
# в сравнение с базовым на вероятность отличного здоровья
data.bachelor <- data
data.bachelor$bachelor <- 0
data.bachelor$master <- 0
ame.bachelor.1 <- mean(predict(model1, group = 4, type = "prob", 
                               me = "bachelor", eps = c(0, 1),
                               newdata = data.bachelor))
ame.bachelor.2 <- mean(predict(model2, group = 4, type = "prob", 
                               me = "bachelor", eps = c(0, 1),
                               newdata = data.bachelor))
ame.bachelor.3 <- mean(predict(model3, group = 4, type = "prob", 
                               me = "bachelor", eps = c(0, 1),
                               newdata = data.bachelor))
c(probit = ame.bachelor.1, logit = ame.bachelor.2, 
  semiparametric = ame.bachelor.3)

# Посчитаем средний предельный эффект наличия магистерского образования
# в сравнение с бакалаврским на вероятность отличного здоровья
ame.mb.1 <- ame.master.1 - ame.bachelor.1
ame.mb.2 <- ame.master.2 - ame.bachelor.2
ame.mb.3 <- ame.master.3 - ame.bachelor.3
c(probit = ame.mb.1, logit = ame.mb.2, semiparametric = ame.mb.3)

# Протестируем гипотезу о том, что степень бакалавра и магистра оказывают
# одинаковый эффект на здоровье
fn_test <- function(object)
{
  coef_val <- coef(object, eq = 1)
  val <- coef_val["master"] - coef_val["bachelor"]
  return(val)
}
test.1 <- delta_method(model1, fn = fn_test)
summary(test.1)
test.2 <- delta_method(model2, fn = fn_test)
summary(test.2)
test.3 <- delta_method(model3, fn = fn_test)
summary(test.3)