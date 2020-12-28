#Подключаем библиотеки
library(readxl)
library(tseries)
library(forecast)
library(ggplot2)
library(aTSA)
library(fGarch)
library(plotly)
#Определяем функции
##Функция полиномиального сглаживания
poly_smooth = function(series, m, p){ 
  T = NROW(series)
  sa = array(NA,dim = T)
  for (i in (m + 1):(T - m)) {
    t = (-m):m
    temp = series[(i - m):(i + m)]
    fit = lm(temp~poly(t,p,raw = TRUE))
    sa[i] = as.numeric(fit$coefficients[1])
    if (i == m + 1) {
      sa[1:m] = predict(fit,data.frame(t = (-m):(-1)))
    }
    if (i == (T - m)) {
      sa[(T - m + 1):T] = predict(fit,data.frame(t = (1):m))
    }
  }
  sa
}
##Функция экспоненциального сглаживания
exp_smooth = function(series,alpha) {
  if (alpha > 1 | alpha < 0) {
    print('Беды с башкой?')
  } else {
    T = NROW(series)
    s = series[1]
    for (i in 1:T) {
      s = c(s,alpha*series[i] + (1 - alpha)*s[i])
    }
    s = s[-1]
    s
  }
}
##Функция выдает датафрейм с альфами, бетами и интенсивностями для гармоник
fourier = function(series) {
  N = NROW(series)
  j = 0:(N %/% 2)
  f = 2*pi*(1:N)
  alpha = c(sum(series)/N)
  intensive = c((alpha[1]**2)*N/2)
  betta = c(0)
  for (i in 1:(N %/% 2)) {
    alpha = c(alpha,2/N*sum(series*cos(f*i/N)))
    betta = c(betta,2/N*sum(series*sin(f*i/N)))
    intensive = c(intensive,(alpha[i + 1]**2 + betta[i + 1]**2)*N/2)
  }
  if (N %% 2 == 0) {
    alpha = c(alpha,sum((-1)**(1:N)*series)/N)
    intensive = c(intensive,alpha[N %/% 2 + 1]**2)
    betta = c(betta,0)
    j = c(j,N %/%  2+1)
  }
  data.frame(alpha,betta,intensive,j)
}
##Функция строит ряд из q основных гармоник на основе результата прошло функции
fourier2 = function(q,spect){
  if(q > 0) {
    spect = spect[order(spect$intensive,decreasing = TRUE),]
    x = c()
    for (i in 1:T) {
      x[i] = sum(spect$alpha[1:(q)]*cos(2*pi*i*spect$j[1:(q)]/T) + spect$betta[1:(q)]*sin(2*pi*i*spect$j[1:(q)]/T))
    }
  }
  x
}
#Загружаем дынные
TimeSeries = read_excel('/Users/nikita_tililitsin/Desktop/U.S. Air Carrier Traffic Statistics - Revenue Passenger Miles.xls', sheet = 'Лист1')
names(TimeSeries) = c('y','date')
TimeSeries$date=as.Date(TimeSeries$date)
TimeSeries$months = months(TimeSeries$date)
T = NROW(TimeSeries)
TimeSeries$t = 1:T
TimeSeries$lny = log(TimeSeries$y)
#Визуализация данных
ggplot(TimeSeries,aes(x = date, y = y)) + 
  geom_line(color="#903749", size = 1)+
theme_minimal()+ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D')) 
ggplot(TimeSeries,aes(x = date, y = lny)) + 
  geom_line(color="#903749", size = 1)+theme_minimal()+ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
#Выявляем тренд
m0 = lm(lny~t,data = TimeSeries)
summary(m0)
ggplot(TimeSeries,aes(x = date, y = lny)) + 
  geom_line(color="#903749", size = 1)+geom_line(aes(y=m0$fitted.values),color="#3F72AF", size = 1)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

m1 = lm(lny~t+sqrt(t),data = TimeSeries)
summary(m1)
ggplot(TimeSeries,aes(x = date, y = lny)) + 
  geom_line(color="#903749", size = 1)+geom_line(aes(y=m1$fitted.values),color="#3F72AF", size = 1)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

m2 = lm(lny~poly(t,2, raw = TRUE),data = TimeSeries)
summary(m2)
ggplot(TimeSeries,aes(x = date, y = lny)) + 
  geom_line(color="#903749", size = 1)+geom_line(aes(y=m2$fitted.values),color="#3F72AF", size = 1.)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

#Выявляем сезонность при помощи фиктивных переменных
ms = lm(lny~t+sqrt(t)+factor(months),data = TimeSeries)
  summary(ms)
ggplot(TimeSeries,aes(x = date, y = lny)) + 
  geom_line(color="#903749", size = 1.)+geom_line(aes(y=ms$fitted.values),color="#3F72AF", size = 1.)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
  
#График остатков и тест Дики-Фулера (не знаю зачем, но пусть будет)
ggplot(TimeSeries,aes(x = date, y = ms$residuals)) + 
  geom_line(color="#903749", size = 1)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
adf.test(ms$residuals)

#Прогноз по тренду и сезонности
newmonths = seq(as.Date(TimeSeries$date[T]),by = 'month', length.out = 13)
newdata = data.frame(t = (T+1):(T+12), months = months(newmonths[-1]))
pr = data.frame(predict(ms,newdata, interval = "prediction", level = 0.95))
fit.val = c(ms$fitted.values,pr$fit)
ggplot() + 
  geom_line(aes(x =TimeSeries$date, y = TimeSeries$lny),color="#903749", size = 1.)+
  geom_line(aes(x=c(TimeSeries$date,newmonths[-1]),y=fit.val),color="#3F72AF", size = 1.)+
  geom_line(aes(x = newmonths[-1],y=pr$upr),color="#FF9A00", size = 1.)+
  geom_line(aes(x = newmonths[-1],y=pr$lwr),color="#FF9A00", size = 1.)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

# ACF and PACF
acf = Acf(ms$residuals)
ggplot() + 
  geom_line(aes(x = acf$lag[,,1], y = acf$acf[,,1]),color="#3F72AF", size = 1)+
  geom_hline(aes(yintercept=0),size = 1,color = "#903749")+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
pacf = pacf(ms$residuals)
ggplot() + 
  geom_line(aes(x = pacf$lag[,,1], y = pacf$acf[,,1]),color="#3F72AF", size = 1)+
  geom_hline(aes(yintercept=0),size = 1,color = "#903749")+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
#Критерий спирмена
spear = cor(ms$residuals,1:T, method = 'spearman')
spear
spear*sqrt((T-2)/(1-spear**2))

#полиномиальноe сглаживаниe
# m = 4 p = 1
sa1 = poly_smooth(TimeSeries$y,4,1)
# m = 4 p = 2
sa2 = poly_smooth(TimeSeries$y,4,2)
ggplot(TimeSeries,aes(x=date,y = y)) + 
  geom_line(color="#903749", size = 1)+
  geom_line(aes(y = sa1),color="#68006C", size = 1)+
  geom_line(aes(y=sa2),color="#00B8A9", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
  

#экспоненциальноe сглаживаниe
s = exp_smooth(TimeSeries$y,.3)
ggplot(TimeSeries,aes(x=date,y = y)) + 
  geom_line(color="#903749", size = 1)+
  geom_line(aes(y = s),color="#00b8a9", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

# Остатки после удаления тренда(т.е. с сезонностью)
# дальше будем моделировать сезонность с помощью гармонического анализа
ggplot(TimeSeries,aes(x=date,y = m1$residuals)) + 
  geom_line(color="#903749", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))


#Гармонический анализ
spect = fourier(m1$residuals)
ggplot(spect,aes(j,intensive))+geom_col(fill='#903749')+theme_minimal()


x = fourier2(4,spect)
ggplot(TimeSeries,aes(x=date,y = m1$residuals)) + 
  geom_line(color="#903749", size = 1)+
  geom_line(aes(y = x),color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
  
ggplot(TimeSeries,aes(x=date,y = lny)) + 
  geom_line(color="#903749", size = 1)+
  geom_line(aes(y = x+m1$fitted.values),color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

#АРИМА
fit = auto.arima(ms$residuals[1:(T-10)])
summary(fit)
ggplot() + 
  geom_line(aes(x=TimeSeries$date[1:(T-10)],y = ms$residuals[1:(T-10)]),color="#903749", size = 1)+
  geom_line(aes(x=TimeSeries$date[1:(T-10)],y = fit$fitted),color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

adf.test(fit$residuals)

arma1 = Arima(ms$residuals[1:(T-10)], order = c(0,0,1))
ggplot() + 
  geom_line(aes(x=TimeSeries$date[1:(T-10)],y = ms$residuals[1:(T-10)]),color="#903749", size = 1)+
  geom_line(aes(x=TimeSeries$date[1:(T-10)],y = arma1$fitted),color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
summary(fit)
summary(arma1)
prar1 = forecast::forecast(arma1,10)
prar2 = forecast::forecast(fit,10)
autoplot(prar1)+geom_line(aes(x=1:T,y=ms$residuals),color='#903749',size=0.5)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
autoplot(prar2)+geom_line(aes(x=1:T,y=ms$residuals),color='#903749',size=0.5)+theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))

cor(TimeSeries$lny[1:(T-10)], ms$fitted.values[1:(T-10)] + fit$fitted)^2

ggplot(TimeSeries,aes(x=date,y=lny)) +
  geom_line(color="#903749",size=1)+
  geom_line(aes(y = ms$fitted.values+c(fit$fitted,prar2$mean)),color="#3F72AF", size = 1.)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))


#Гарч
trend=lm(lny~poly(t,2,raw = TRUE), data=TimeSeries[1:(T-10),])
mgarch=garchFit(formula = ~arma(2,0)+garch(1,1),data = trend$residuals)
summary(mgarch)
predgarch=predict(mgarch,10, conf = .95, plot=TRUE)
fullfit=trend$fitted.values+mgarch@fitted
plot(TimeSeries$lny, type='l')
lines(trend$fitted.values)
lines(fullfit, col = 'red')
predtrend=predict.lm(object = trend,newdata = data.frame(t =(T-9):T),interval = "prediction",level = 0.95)
f=plot_ly(
  y=TimeSeries$lny,
  x=TimeSeries$t,
  type='scatter',
  mode='lines',
  name = 'Исходный ряд'
) %>%
  layout(
    title = 'Тренд + ARMA + Garch',
    xaxis = list(
      title = 't'
    ),
    yaxis = list(
      title = 'Value'
    )
  ) %>%
  add_trace(y=trend$fitted.values, x=1:(T-10), name = 'Тренд') %>%
  add_trace(y=fullfit, x=1:(T-10), name = 'Тренд+ARMA+Garch') %>%
  add_trace(y=predgarch$meanForecast + predtrend[,1],x=(T-9):T, name = 'Прогноз')
f

#Коинтеграция
Coint = read_excel('/Users/nikita_tililitsin/Desktop/U.S. Air Carrier Traffic Statistics - Revenue Passenger Miles.xls', sheet = 'Лист 2')
Coint$t = 1:NROW(Coint)
Coint$RUB=as.numeric(Coint$RUB)
Coint$KZT=as.numeric(Coint$KZT)
ggplot(Coint,aes(x = Date, y = RUB)) + 
  geom_line(color="#903749", size = 1)+
  theme_minimal()
ggplot(Coint,aes(x = Date, y = KZT)) + 
  geom_line(color="#3F72AF", size = 1)+
  theme_minimal()

ndiffs(Coint$RUB)
ndiffs(Coint$KZT)
KZT_1 = diff(Coint$KZT,differences = 1)
RUB_1 = diff(Coint$RUB,differences = 1)
ggplot(Coint[-1,],aes(x = Date, y = RUB_1)) + 
  geom_line(color="#903749", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
ggplot(Coint[-1,],aes(x = Date, y = KZT_1)) + 
  geom_line(color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
tseries::adf.test(KZT_1)
tseries::adf.test(RUB_1)



mcoint=lm(RUB~KZT-1,data = Coint)
summary(mcoint)

mcoint1 = lm(mcoint$residuals[-1]~mcoint$residuals[-131]-1)
summary(mcoint1)

coint.test(Coint$KZT, Coint$RUB,d=1)

ecm = lm(RUB_1~KZT_1+mcoint$residuals[-131]-1)
summary(ecm)
  ggplot(Coint[-1,],aes(x = Date, y = RUB_1)) + 
  geom_line(color="#903749", size = 1)+
  geom_line(aes(y=ecm$fitted.values),color="#3F72AF", size = 1)+
  theme_minimal()+
  ylab("") +xlab("")+theme(text = element_text(size=15),axis.text.x = element_text(color= '#48466D'),axis.text.y = element_text(color= '#48466D'))
           
