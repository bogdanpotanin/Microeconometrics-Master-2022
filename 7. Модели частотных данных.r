# --------
# Потанин Богдан Станиславович
# Микроэконометрика в R :)
# Семинар 7. Модель с частотными данными
# --------

# Отключим scientific notation
options(scipen = 999)

# Подключим дополнительные библиотеки
library("MASS")                                            # отрицательная биномиальная
                                                           # регрессия
library("pscl")                                            # модель с инфляцией нулей

# Будем моделировать посещение индивидом врача
set.seed(123)
n <- 10000
income <- runif(n)
age <- runif(n)
beta <- c(0.1, 0.8, 1)
lambda <- exp(beta[1] + beta[2] * income + beta[3] * age)

# Пуассоновская регрессия
doctor <- rpois(n, lambda)
model.pois <- glm(doctor ~ income + age, family = "poisson")
summary(model.pois)
# предсказания
head(predict(model.pois))

# Отрицательная биномиальная регрессия
theta <- 0.3
eps <- log(rgamma(n, shape = 1 / theta, scale = theta))
lambda <- exp(beta[1] + beta[2] * income + beta[3] * age + eps)
doctor <- rpois(n, lambda)
model.nb <- glm.nb(doctor ~ income + age)
summary(model.nb)
# предсказания
head(predict(model.nb))
# Вместо theta пакет оценивает 1 / theta
theta_hat <- 1 / model.nb$theta

# Модель с инфляцией нулей
gamma <- c(0.1, 3, -2)
u <- rnorm(n)
z.li <- gamma[1] + gamma[2] * income + gamma[3] * age 
z.star <- z.li + u
z <- rep(1, n)
z[z.star < 0] <- 0
doctor[z == 0] <- 0
# но в пакете оценивает, меняя местами 0 и 1, поэтому
# в уравнении нулей коэффициенты будут с обратными знаками
model.zim <- zeroinfl(doctor ~ income + age | 
                               income + age, 
                      dist = "negbin", link = "probit")
summary(model.zim)
# вместо теты вновь оценивается 1 / theta
theta_hat <- 1 / model.zim$theta
# предсказания
head(predict(model.zim ))
