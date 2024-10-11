# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 9. Панельные модели бинарного выбора
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("numDeriv")                                      # численное дифференцирование
library("GJRM")                                          # оценивание систем
                                                         # бинарных уравнений
library("mnorm")                                         # двумерное нормальное распределение
library("hpa")                                           # полу-непараметрический подход
library("ggplot2")                                       # красивые графики
library("np")                                            # ядерное оценивание параметров
                                                         # моделей бинарного выбора
library("Ecdat")                                         # встроенные данные
library("margins")                                       # расчет предельных эффектов
library("pglm")                                          # панельные модели бинарного выбора
library("sandwich")                                      # кластеризация
library("lmtest")                                        # тестирование гипотез

#---------------------------------------------------
# Часть 1. Генерация данных
#---------------------------------------------------

# Симулируем данные
set.seed(123)                                            # для воспроизводимости
n       <- 1000                                          # число наблюдений    
m       <- 10                                            # число наблюдений по каждой категории
X       <- rmnorm(n,                                     # симулируем n наблюдений из многомерного
                                                         # нормального распределения
                  c(0, 0, 0),                            # с нулевым вектором математических ожиданий и
                  matrix(c(1,    0.2,  0.3,              # следующей ковариационной матрице
                           0.2,    1, -0.1,
                           0.3, -0.1,    1),
                   ncol = 3, byrow = FALSE))
X[, 2]  <- X[, 2] > 0.5                                  # сделаем регрессор X2 бинарным
X       <- cbind(1, X)                                   # добавим константу как дополнительный регрессор
mu      <- 0                                             # параметры распределения
sigma   <- 1                                             # случайных ошибок
u       <- rnorm(n,                                      # случайные ошибки из нормального распределения
           mu,                                           # с математическим ожиданием mu и
           sigma)                                        # стандартным отклонением sigma
mu_a    <- 0
sigma_a <- 0.8
a       <- rnorm(n / m, mean = mu_a, sd = sigma_a)       # симулируем случайные эффекты
a       <- rep(a, each = m)                              # распределяем случайные эффекты
                                                         # по наблюдениям

# Соберем регрессоры в датафрейм
h <- data.frame("id" = rep(1:(n / m), each = m),         # уникальный идентификатор наблюдения
                "Intercept" = X[, 1],                    # константа
                "skills" = X[, 2],                       # показатель, отражающий уровень
                # развития профессиональных равыков
                "male" = X[, 3],                         # пол: мужчина = 1, женщина = 0
                "experience" = X[, 4],                   # показатель, отражающий опыт работы
                "weight" = runif(n))                     # показатель веса, симулированный из
                                                         # стандартного равномерного распределения

# Создадим линейный предиктор, который определяет
# зависимость латентной переменной от различных
# независимых переменных
gamma <- c(0.5, 0.8, 0.7, 0.8, -0.3, 0.5, 0)             # оцениваемые регрессионные коэффициенты
z_li  <- gamma[1] * h$Intercept +                        # константа
         gamma[2] * h$skills +                           # навыкы
         gamma[3] * h$male +                             # брак
         gamma[4] * h$experience +                       # опыт
         gamma[5] * h$experience ^ 2 +                   # квадрат опыта
         gamma[6] * h$skills * h$male +                  # взаимодействие навыков и пола
         gamma[7] * h$weight                             # показатель веса
z_star <- z_li + u + a                                   # латентная переменная

# Создадим наблюдаемую зависимую переменную,
# отражающую работает индивид или нет
z      <- as.numeric(z_star >= 0)                        # наблюдаемое значение переменной
z      <- matrix(z, ncol = 1)                            # как матрица с одним столбцом 
h$work <- z                                              # добавим в датафрейм переменную 
                                                         # на трудоустройство
sum(h$work  / n)                                         # доля работающих индивидов

# Посмотрим на данные
head(h)

#---------------------------------------------------
# Часть 2. Оценивание параметров модели
#---------------------------------------------------

# Для начала воспользуемся обычной пробит моделью
model_probit <- glm(formula = work ~ skills + male + 
                                     experience + 
                                     I(experience ^ 2) +  
                                     I(skills * male) +
                                     weight,                                                 
                    data = h,
                    family = binomial(link = "probit"))      
gamma_probit <- coef(model_probit)                       # сохраним коэффициенты
summary(model_probit)                                    # посмотрим на результат

# Применим модель со случайными эффектами
model_panel <- pglm(work ~ skills +                      # формула по 
                           male +                        # аналогии с glm
                           experience +
                           I(experience ^ 2) +
                           I(skills * male) +
                           weight,
                    index = "id",                        # группирующая переменная
                    model = "random",                    # тип эффекта
                    method = "bfgs",                     # метод оптимизации
                    R = 20,                              # чем, больше, тем точней,
                                                         # но дольше расчеты, за счет
                                                         # использования большего числа весов
                                                         # в разложении квадратуры Гаусса - Эрмита.
                    data = h,                            # данные по аналогии  glm
                    
                    family = binomial(link = 'probit'))  # тип модели по аналогии с glm

# Посмотрим на результат
summary(model_panel)                                     # общая выдача
sigma_a_est <- model_panel$estimate["sigma"] / sqrt(2)   # оценка стандартного отклонения
                                                         # случайного эффекта
gamma_panel <- coef(model_panel)                         # оценки коэффициентов из которых
gamma_panel <- gamma_panel[names(gamma_panel) !=         # для удобства исключим sigma_a
                           "sigma"]

# Сравним точность оценок
sd_panel <- sqrt(sigma ^ 2 + sigma_a ^ 2)
coef_df <- cbind("Real"             = gamma,
                 "Probit"           = gamma_probit,
                 "Probit Adj"       = gamma_probit * sd_panel,
                 "Panel"            = gamma_panel,
                 "Probit Error"     = (gamma - gamma_probit) ^ 2,
                 "Probit Adj Error" = (gamma - gamma_probit * sd_panel) ^ 2,
                 "Panel Error"      = (gamma - gamma_panel) ^ 2)
print(round(coef_df, 4))
# Важно - при сравнение оценок коэффициентов обычной и панельной пробит модели
#         важно учитывать различия в дисперсиях случайных ошибок

# Сравним модели по информационным критериям
data.frame("Probit" = AIC(model_probit),
           "Panel"  = AIC(model_panel))

#---------------------------------------------------
# Часть 3. Оценивание условных вероятностей
#---------------------------------------------------

# Оценим вероятность занятости для Бориса
Boris <- data.frame("work" = 1,
                    "skills" = 0.2,                      # укажем характеристики
                    "male" = 1,                          # Бориса в датафрейме
                    "experience" = -0.3,
                    "weight" = 0.2)
  # получим матрицу характеристик Бориса в
  # соответствии с формулой, по которой
  # рассчитывается линейный предиктор
X_Boris <- model.frame(formula = formula(model_panel), 
                       data = Boris)[-1]
  # рассчитаем вероятность по обычной пробит модели
prob_Boris_probit <- predict(model_probit, 
                             newdata = Boris,
                             type = "response")
  # рассчитаем вероятность по модели со 
  # случанйыми эффектами
li_Boris <- gamma_panel[1] +                             # линейный предиктор
            sum(X_Boris * gamma_panel[-1])               
prob_Boris_panel <- pnorm(li_Boris,                      # вероятность
                          sd = sqrt(1 + sigma_a_est ^ 2))
  # сравним оценки вероятностей занятости
data.frame("Probit" = prob_Boris_probit,
           "Panel" = prob_Boris_panel)

# Рассчитаем вероятность занятости Бориса,
# используя закон больших чисел ЗБЧ:
# P(z = 1) = P(e + a > - xb) = P(e > -a - xb) = 
# E(P(e > -a - xb)) = E[E(P(e > -a - xb) | a)] =
# E[E(F(xb + a) | a)], 
# где E(F(xb + a) | a), при не фиксированном "a",
# это случайная величина, которая зависит от
# случайной величины "a", что позволяет
# применить ЗБЧ
n_sim <- 1000000                                         # чем больше, чем точней
a_new <- rnorm(n_sim, sd = sigma_a_est)                  # выборка из того же распределения,
                                                         # что и у случайного эффекта
prob_Boris_LLN <- mean(pnorm(li_Boris + a_new))          # оцениваем вероятность
  # сравним результат
data.frame("Analytical " = prob_Boris_panel,
           "LLN"         = prob_Boris_LLN)

# Теперь рассмотрим Бориса во втором
# периоде времени
Boris_new <- data.frame("work" = 1,
                        "skills" = 0.3,
                        "male" = 1,
                        "experience" = -0.1,
                        "weight" = 0.5)
X_Boris_new <- model.frame(formula = formula(model_panel), 
                           data = Boris_new)[-1]
li_Boris_new <- gamma_panel[1] + sum(X_Boris * gamma_panel[-1])

# Используя ЗБЧ оценим вероятность того, что в
# Борис был занят в обоих периодах времени:
# P(z1 = 1, z2 = 1) = P(e1 > -a - x1b, e2 > -a - x2b) = 
# E(P(e1 > -a - x1b, e2 > -a - x2b))
# E[E(P(e1 > -a - x1b, e2 > -a - x2b) | a)] =
# E[E(P(e1 > -a - x1b) * P(e2 > -a - x2b) | a)]
p_11_LLN <- mean(pnorm(li_Boris + a_new) *
                 pnorm(li_Boris_new + a_new))

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 3.1.    Изучите, как зависит преимущество в точности оценок
#         пробит модели со случайными эффектами от дисперсии
#         случайных эффектов: сравните результаты при
#         дисперсиях 0.1, 0.5, 1, 2, и 5.
# 3.2.    Оцените вероятность занятости для индивида
#         с произвольными характеристиками
# 3.3.    Оцените предельный эффект:
#         1*)   опыта работы на вероятность занятости у Бориса
#         2*)   стандартного отклонения случайного эффекта
#               на вероятность занятости у Бориса
#         3**)  опыта работы на вероятность того, что в первый
#               период времени Борис работал, а во воторой - нет
#         4***) опыта работы на вероятность того, что Борис
#               работает во второй периодб при услови, что он
#               не работал в первый период
# 3.4**.  Проверьте гипотезу о равенстве нулю
#         предельного эффекта опыта работы
#         на вероятность занятости у Бориса
# 3.5.    Постройте 90%-й доверительный интервал для:
#         1**)  вероятности занятости Бориса
#         2***) вероятности занятости Бориса два периода
#               времени подряд

#---------------------------------------------------
# Часть 4. Кластеризация стандартных ошибок
#---------------------------------------------------

# В качестве альтернативы модели со случайными эффектами
# можно использовать обычную пробит модель, оценки которой
# состоятельны (с учетом корректировки на дисперсию), 
# но требуют учета панельной структуры для расчета стандартных ошибок

# Оценим асимптотическоую ковариационную матрицу с учетом кластеризации
vcov_cluster <- vcovCL(model_probit, sandwich = TRUE, type = "HC0", 
                       cadjust = FALSE, cluster = ~ id)

# Протестируем гипотезы с оцененной ковариационной матрицей
coeftest(model_probit, vcov = vcov_cluster)

# Сравним стандартные ошибки
se_df <- cbind("Probit no Cluster" = sqrt(diag(vcov(model_probit))),
               "Probit Cluster"    = sqrt(diag(vcov_cluster)))
print(round(se_df, 4))

# Начнет работать в версии 2.0.1 пакета switchSelection
# model_probit2 <- msel(formula(model_probit), data = h, cluster = ~ id)
# summary(model_probit2)

#---------------------------------------------------
# Часть 5. Применение панельной пробит модели
#          к реальным данным
#---------------------------------------------------

# По мотивам статьи:
# Vella, F. and M. Verbeek (1998) “Whose wages do unions raise? 
# A dynamic model of unionism and wage”, Journal of Applied 
# Econometrics, 13, 163–183.

# Воспользуемся датафреймом с информацией о зарплатах
# и членстве в профсоюзе
data(UnionWage)
h <- as.data.frame(UnionWage)

# Для удобства перекодируем переменные
h$year <- as.numeric(h$year)
  # переменная на здоровье            
levels(h$health) <- c(0, 1)                              # изменяем значения факторов на 0 и 1
h$health <- as.numeric(levels(h$health))[h$health]       # преобразуем факторную переменную в численную
  # переменная на проживание в сельской местности           
levels(h$rural) <- c(0, 1)                               # изменяем значения факторов на 0 и 1
h$rural <- as.numeric(levels(h$rural))[h$rural]          # преобразуем факторную переменную в численную
  # переменная на членство в профсоюзе       
levels(h$union) <- c(0, 1)                               # изменяем значения факторов на 0 и 1
h$union <- as.numeric(levels(h$union))[h$union]          # преобразуем факторную переменную в численную
  # индивид характеризует себя как белый американец
h$white <- 0                                             # создаем переменную с 0 значениями
h$white[h$com == "white"] <- 1                           # изменяем 0 на 1 для нужных индивидов
  # индивид характеризует себя как афроамериканец
h$black <- 0                                             # создаем переменную с 0 значениями
h$black[h$com == "black"] <- 1                           # изменяем 0 на 1 для нужных индивидов
  # индивид характеризует себя как испанец
h$hisp <- 0                                              # создаем переменную с 0 значениями
h$hisp[h$com == "hispanic"] <- 1                         # изменяем 0 на 1 для нужных индивидов
  # членство в профсоюзе с прошлом году
h$union_lag <- 0
  #h$union_lag[h$year == 1980] <- NA
h$union_lag[h$year >= 1981] <- h$union[h$year < 1987] 

# Оценим вероятность членства в профсоюзе
# используя обычную пробит модель
model_probit <- glm(union ~ I(log(exper + 1)) +
                            school +
                            married +
                            black + hisp +
                            rural +
                            health + 
                            region + sector +
                            as.factor(year),
                   data = h,             
                   family = binomial('probit'))   
summary(model_probit)
coef_probit <- coef(model_probit)

# Оценим вероятность членства в профсоюзе
# используя пробит модель со случайными
# эффектами
model_RE <- pglm(union ~ I(log(exper + 1)) +
                         school +
                         married +
                         black + hisp +
                         rural +
                         health + 
                         region + sector +
                         as.factor(year),
                    index = "id",
                    model = "random",
                    method = "bfgs",
                    data = h,             
                    
                    family = binomial('probit'))                   
summary(model_RE)                                        # общая выдача
sigma_a_est <- model_RE$estimate["sigma"] / sqrt(2)      # оценка стандартного отклонения
                                                         # случайного эффекта
gamma_RE <- coef(model_RE)                               # оценки коэффициентов из которых
gamma_RE <- gamma_RE[names(gamma_RE) != "sigma"]         # для удобства исключим sigma_a

# ЗАДАНИЯ (* - непросто, ** - сложно, *** - брутально)
# 5.1.    Сравните модель со случайными эффектами с
#         моделью с гетероскедастичной случайной
#         ошибкой по:
#         1)     информационным критериям
#         2)     предсказанным вероятностям
#         3*)    предсказанным предельным эффектам
# 5.2.    Сравните пробит и логи модели со случайными
#         эффектами по:
#         1)     информационным критериям
#         2)     предсказанным вероятностям
#         3*)    предсказанным предельным эффектам

#---------------------------------------------------
# Часть 6. Реализация пробит модели со случайными
#          эффектами при помощи симуляционного
#          метода максимального правдоподобия
#---------------------------------------------------

# Логарифм функции правдоподобия, максимизируемый
# при использовании панельной пробит регрессии
ProbitRELnL <- function(x, z, X, id, n_sim = 1000)       # функция правдоподобия
{
  sigma_a <- x[1]                                        # стандартное отклонение
                                                         # случайных эффектов
  
  if(sigma_a <= 0)                                       # избегаем отрицательных
  {                                                      # стандартных отклонений
    return(-99999999999999999)
  }
  
  gamma <- matrix(x[-1], ncol = 1)                       # вектор регрессионных коэффициентов
                                                         # регрессионных коэффициентов,
                                                         # переводим в матрицу с одним столбцом
  z_est <- X %*% gamma                                   # оценка математического ожидания 
                                                         # латентной переменной

  n_obs <- nrow(X)                                       # количество наблюдений
  
  is_z_0 <- z == 0                                       # вектор условий z = 0
  is_z_1 <- z == 1                                       # вектор условий z = 1
  
  a_new <- rnorm(n_sim, sd = sigma_a)                    # симулируем значения из того же
                                                         # распределения, что и у
                                                         # случайных эффектов
  
  lnL <- 0                                               # логарифм функции правдоподобия
         
  P_i <- rep(1, n_sim)                                   # совместная вероятность для
                                                         # индивида из категории

  for(i in unique(id))                                   # для каждой категории
  {
    for(j in which(id == i))                             # для каждого наблюдения
    {                                                    # в данной категории
      z_li_a <- z_est[j] + a_new                         # считаем сумму линейного индекса
                                                         # и случайного эффектоа
      if(is_z_1[j])                                      # если это успех, то считаем
      {                                                  # вероятность как для успеха
        P_i <- P_i * pnorm(z_li_a)                       # симулированных значений
      } else {                                           # если это не успех, то считаем
        P_i <- P_i * (1 - pnorm(z_li_a))                 # при всех реализациях
      }                                                  # симулированных значений
    }
    lnL <- lnL + log(mean(P_i))                          # добавляем посчитанный результат
                                                         # к логарифму функции правдоподобия
    P_i <- rep(1, n_sim)                                 # возвращаем вероятность к 1 перед
                                                         # переходом к наблюдениям из 
                                                         # новой категории
  }
  
  cat(paste("sigma_a = ", sigma_a, "\n"))
  cat(paste("lnL = ", lnL, "\n"))
  
  return(lnL)
}

ProbitRE <- function(formula,                            # формула
                     data,                               # датафрейм содержащий переменную id,
                                                         # по которой группируются данные
                     x0 = NULL,                          # начальная точка
                     n_sim = 1000)                       # число симуляций для ЗБЧ
                                                         
{
  d <- model.frame(formula, data)                        # извлекаем переменные согласно формуле
  id <- data$id                                          # достаем группирующую переменную из данных
  
  z <- as.matrix(d[, 1], ncol = 1)                       # зависимая переменная как первый
                                                         # столбец в d
  X <- as.matrix(d[, -1])                                # независимые переменные как все переменные
                                                         # из data кроме зависимой
  x0_n <- ncol(X) + 1                                    # число оцениваемых параметров
  
  X_names <- names(d)[-1]                                # имена независимых переменных
  
  
  if(is.null(x0))
  {
    x0 <- rep(0, x0_n)
    model_probit <- glm(z~. + 0,                         # берем начальные точки
                        data = as.data.frame(X),         # из обычной пробит регрессии
                        family = binomial('probit'))
    x0[-1] <- coef(model_probit)
    x0[1] <- 0.01
  }

  result <- optim(par = x0,                              # в качестве начальных точек возьмем нули
                  method = "BFGS",                       # численный метод оптимизации
                  fn = ProbitRELnL,                      # максимизируемая функция правдоподобия
                  control = list(maxit = 10000000,       # чтобы минимизационную задачу превратить
                                 fnscale = -1,           # в максимизационную умножаем функцию на -1
                                 reltol = 1e-10),        # установим достаточно высокую точность          
                  hessian = TRUE,                        # вернем Гессиан функции
                  X = X, z = z,                          # аргументы оптимизируемой функции 
                  id = id, n_sim = n_sim)                 

  
  gamma_est <- result$par[-1]                            # оценки коэффициентов
  sigma_a_est <- result$par[1]                           # оценки стандартного отклонения
                                                         # случайного эффекта
  names(gamma_est) <- X_names                            # сопоставляем имена для оценок коэффициентов
  
  as_cov_est <- -solve(result$hessian)                   # оценка асимптотической ковариационной
  colnames(as_cov_est) <- c("sigma_a", X_names)          # матрицы полученных оценок
  rownames(as_cov_est) <- c("sigma_a", X_names)          # сопоставляем имена
  
  return_list <- list("gamma" = gamma_est,               # возвращаем оценки коэффициентов и
                      "cov" = as_cov_est,                # асимптотической ковариационной матрицы
                      "data" = data,                     # возвращаем использовавшийся датафрейм
                      "lnL" = result$value,              # возвращаем логарифм функции правдоподобия
                      "X" = X,                           # возвращаем матрицу регрессоров
                      "z" = z,                           # возвращаем зависимую переменную  
                      "sigma_a" = sigma_a_est,           # оценка стандартного отклонения
                                                         # случайных эффектов
                      "formula" = formula)               # возвращаем использовавшуюся формулу                          
  
  class(return_list) <- "probitRE"                       # для удобства назначим класс
                                                         # возвращаемой из функции переменной
  return(return_list)                                    # возвращаем результат                               
}
# Воспользуемся созданной функцией и применим
# метод, именуемый maximum simulated likelihood
model <- ProbitRE(work ~ Intercept +                     # указываем формулу с константой, поскольку
                  skills + male +                        # она не учитывается автоматически
                  experience + I(experience ^ 2) +                   
                  I(skills * male) +
                  weight,
                 data = h,
                 n_sim = 1000)                                    
gamma_est <- model$gamma                                 # получаем оценки коэффициентов
gamma_cov_est <- model$cov                               # получаем оценку асимптотической

# Проверка за счет
# x0 = c(model_panel$estimate[8] / sqrt(2),
#        model_panel$estimate[1:7])