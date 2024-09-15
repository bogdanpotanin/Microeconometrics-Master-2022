# --------
# Фамилия, Имя, Отчество
# Микроэконометрика качественных данных
# Домашнее задание
# --------

# Отключим scientific notation
options(scipen = 999)

#---------------------------------------------------
# Симуляция данных
#---------------------------------------------------

# Рассмотрим факторы, влияющие на вероятность того, что
# индивиду понравится ваша стратегическая игра про
# кошачий город

library("mnorm")
library("switchSelection")

# Симулируем данные
set.seed(888)                                              # для воспроизводимости
n        <- 5000                                           # число индивидов в выборке
h        <- data.frame(ind = 1:n)
h$income <- round(exp(rnorm(n, 10, 0.7)) / 1000)* 1000     # доход
h$age    <- round(runif(n) ^ 2 * 50 + 18)                  # возраст
h$hours  <- round(rbeta(n, 1, 10) * 100 + 1)               # играет часов в неделю
h$health <- as.factor(sample(c("good", "medium", "bad"),   # здоровье
                             n, replace = TRUE, 
                             prob = c(0.3, 0.5, 0.2)))  
h$marriage <- rbinom(n, 1, 0.7)                            # состоит в браке
h$chl      <- rbinom(n, 1, 0.7)                            # факт наличия детей
h$male     <- rbinom(n, 1, 0.5)                            # половая принадлежность
h$residence                   <- round(runif(n, 0, 2))     # место проживания
h$residence[h$residence == 0] <- "capital"
h$residence[h$residence == 1] <- "city"
h$residence[h$residence == 2] <- "village"
h$residence                   <- as.factor(h$residence)
h$cat                         <- rbinom(n, 1, 0.6)         # факт наличия кота
h$dog                         <- rbinom(n, 1, 0.5)         # факт наличия собаки
h$bugs                        <- rpois(n, 5)               # число обнаруженных багов
h$price  <- round(rbeta(n, 1, 3) * 200 + 32)               # стоимость компьютера (в тысячх рублей)
h$design <- rbinom(n, 1, 0.5)                              # дизайн кошек и собак в меню выбора рассы
h$educ   <- as.factor(sample(c("basic", "vocational",      # уровень образования
                               "higher", "phd"),     
                             n, replace = TRUE, 
                             prob = c(0.4, 0.3, 0.2, 0.1)))


# Симулируем случайные ошибки из двумерного
# нормального распредеделения
rho <- 0.7                                           # корреляция между
                                                     # случайными ошибками уравнений
eps <- rmnorm(n, 
              mean = c(0, 0), 
              sigma = matrix(c(1, rho, rho, 1), 
                             ncol = 2))
head(eps, 5)
eps1 <- eps[, 1]                                     # случайные ошибки уравнения дефолта
eps2 <- eps[, 2]                                     # случайные ошибки уравнения
                                                     # стабильной работы

# Латентная переменная уравнения удовлетворенности
game_lp <- 0.5 * scale(log(h$income))            +   
             1 * scale(h$age)                    +
          -0.5 * scale(h$age ^ 2)                +
             1 * scale(h$hours)                  +
          -0.6 * scale(h$hours ^ 2)              +
             0 * scale(h$health == "medium")     +
             0 * scale(h$health == "good")       +
           0.4 * scale(h$marriage)               + 
           0.3 * scale(h$chl)                    + 
          -0.4 * scale(h$male)                   + 
           0.3 * scale(h$residence == "capital") + 
           0.5 * scale(h$residence == "village") + 
           0.7 * scale(h$cat)                    + 
           0.6 * scale(h$dog)                    + 
          -0.2 * scale(h$cat * h$dog)            +
          -0.5 * scale(h$bugs)                   + 
           0.6 * scale(log(h$price))             + 
             0 * scale(h$design)                 + 
           0.2 * scale(h$educ == "vocational")   + 
           0.4 * scale(h$educ == "higher")       + 
           0.6 * scale(h$educ == "phd")
game_lp   <- scale(game_lp)
game_star <- game_lp + eps1

# Уровень удовлетворенности
cuts                                                  <- c(-1, 0)
game                                                  <- rep(0, n)
game[(game_star >= cuts[1]) & (game_star <= cuts[2])] <- 1
game[game_star  > cuts[2]]                            <- 2
h$game                                                <- game
table(game)

# Латентная переменная уравнения выбора
choice_lp <- 0.4 * scale(log(h$income))          +
               1 * scale(h$age)                    +
            -0.6 * scale(h$age ^ 2)                +
               1 * scale(h$hours)                  +
               1 * scale(h$hours ^ 2)              +
            -0.1 * scale(h$health == "medium")     +
            -0.2 * scale(h$health == "good")       +
             0.2 * scale(h$marriage)               + 
             0.3 * scale(h$chl)                    + 
            -0.4 * scale(h$male)                   + 
             0.2 * scale(h$residence == "capital") + 
             0.4 * scale(h$residence == "village") + 
               1 * scale(h$cat)                    + 
               1 * scale(h$dog)                    + 
            -0.3 * scale(h$cat * h$dog)            +
            -0.5 * scale(h$bugs)                   + 
             0.6 * scale(log(h$price))             + 
               1 * scale(h$design)                 + 
             0.2 * scale(h$educ == "vocational")   + 
             0.4 * scale(h$educ == "higher")       + 
             0.6 * scale(h$educ == "phd")
choice_lp   <- scale(choice_lp)
choice_star <- choice_lp + eps2

# Выбор
choice   <- as.numeric(choice_star >= -0.3)
h$choice <- choice
table(choice)

# Анализ
model <- msel(formula = list(game ~ I(log(income)) +
                                    age + 
                                    I(age ^ 2 / 100) + 
                                    hours + 
                                    I(hours ^ 2 / 100) + 
                                    I(health == "medium") + 
                                    I(health == "good") + 
                                    marriage + 
                                    chl + 
                                    male + 
                                    I(residence == "capital") + 
                                    I(residence == "city") +
                                    cat + 
                                    dog + 
                                    I(cat * dog) +
                                    bugs + 
                                    I(log(price)) +
                                    I(educ == "vocational") + 
                                    I(educ == "higher") + 
                                    I(educ == "phd") +
                                    choice,
                              choice ~ I(log(income)) +
                                       age + 
                                       I(age ^ 2 / 100) + 
                                       hours + 
                                       I(hours ^ 2 / 100) + 
                                       I(health == "medium") + 
                                       I(health == "good") + 
                                       marriage + 
                                       chl + 
                                       male + 
                                       I(residence == "capital") + 
                                       I(residence == "city") +
                                       cat + 
                                       dog + 
                                       I(cat * dog) +
                                       bugs + 
                                       I(log(price)) +
                                       I(educ == "vocational") + 
                                       I(educ == "higher") + 
                                       I(educ == "phd") +
                                       design), data = h)
summary(model)

#write.xlsx(h, "G:\\Microeconometrics-Tutorial\\Домашнее задание\\Data.xlsx")
#read.xlsx("G:\\Microeconometrics-Tutorial\\Домашнее задание\\Data.xlsx", sheetIndex = 1)
saveRDS(h, file = "E:\\Преподавание\\Микроэконометрика магистратура\\ДЗ БАК 2024-2025\\homework.rds")
