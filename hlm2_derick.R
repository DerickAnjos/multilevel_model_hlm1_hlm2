################################################################################
#               INSTALAÇÃO E CARREGAMENTO DE PACOTES NECESSÁRIOS               #
################################################################################
# Pacotes utilizados (durante os exemplos irei citar os pacotes utilizados)
pacotes <- c("plotly","tidyverse","reshape2","knitr","kableExtra","rgl","car",
             "nlme","lmtest","fastDummies","msm","lmeInfo","jtools")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

# Algoritmo para determinação dos erros-padrão das variâncias no componente de
# efeitos. Fonte: MBA Data Science and Analytics - USP/Esalq

stderr_nlme <- function(model){
  if(base::class(model) != "lme"){
    base::message("Use a lme object model from nlme package")
    stop()}
  resume <- base::summary(model)
  if(base::length(base::names(model$groups))==1){
    m.type <- "HLM2"
  } else if(base::length(base::names(model$groups))==2){
    m.type <- "HLM3"
  }
  if(m.type == "HLM2"){
    vcov_matrix <- model$apVar
    logs_sd_re <- base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re)==2){
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE Components`=base::c("Var(v0j)","Var(e)"),
                                  `Variance Estimatives`= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                                  base::exp(logs_sd_re[[2]])^2),
                                  `Std Err.`=base::c(stderr_tau00,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
    else{
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau01 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(v0j)","Var(v1j)","Var(e)"),
                                  Estimatives= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                       base::exp(logs_sd_re[[2]])^2,
                                                       base::exp(logs_sd_re[[4]])^2),
                                  Std_Err=base::c(stderr_tau00,
                                                  stderr_tau01,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                            base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                                                               base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
  if(m.type == "HLM3"){
    vcov_matrix <- model$apVar
    logs_sd_re <-  base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re) == 3){
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x3)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(t00k)","Var(v0jk)","Var(e)"),
                                  Estimatives=base::c(base::exp(logs_sd_re)[[2]]^2,
                                                      base::exp(logs_sd_re)[[1]]^2,
                                                      base::exp(logs_sd_re)[[3]]^2),
                                  Std_Err=base::c(stderr_tau_u000,
                                                  stderr_tau_r000,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    } 
    else{
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_r100 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u100 <- msm::deltamethod(~exp(x5)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x7)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE_Components`=base::c("Var(t00k)","Var(t10k)",
                                                          "Var(v0jk)","Var(v1jk)",
                                                          "Var(e)"),
                                  `Variance Estimatives`=base::c(base::exp(logs_sd_re)[[4]]^2,
                                                                 base::exp(logs_sd_re)[[5]]^2,
                                                                 base::exp(logs_sd_re)[[1]]^2,
                                                                 base::exp(logs_sd_re)[[2]]^2,
                                                                 base::exp(logs_sd_re)[[7]]^2),
                                  `Std Err.`=base::c(stderr_tau_u000,
                                                     stderr_tau_u100,
                                                     stderr_tau_r000,
                                                     stderr_tau_r100,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                            base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                                                               base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
}

# Carregando a base de dados ----------------------------------------------

load(file = 'estudante_escola.RData')

estudante_escola %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, 
                font_size = 12)

summary(estudante_escola)
glimpse(estudante_escola)

estudante_escola %>% 
  group_by(escola) %>% 
  summarise(quantidade = n()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

estudante_escola %>% 
  group_by(escola) %>% 
  summarise(desempenho_medio = mean(desempenho, na.rm = T)) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

ggplotly(
  ggplot(data = estudante_escola)+
    geom_density(aes(x = desempenho, fill = escola, color = escola), 
                 position = 'identity',alpha = .3, size = .5, )+
    scale_color_viridis_d()+
    scale_fill_viridis_d()+
    theme_classic()
)

estudante_escola %>% 
  group_by(escola) %>% 
  mutate(linhas = 1:n()) %>% 
  mutate(x = unlist(density(desempenho, n = max(linhas))['x']),
         y = unlist(density(desempenho, n = max(linhas))['y'])) %>% 
  ggplot()+
  geom_area(aes(x = x , y = y, group = escola, fill = escola), color = 'black',
            alpha = .3)+
  geom_histogram(aes(x = desempenho, y = ..density.., fill = escola),
                color = 'black', position = 'identity', alpha = .1) + 
  facet_wrap(~escola)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_bw()

ggplotly(
  estudante_escola %>% 
    ggplot(aes(x = horas, y = desempenho))+
    geom_point()+
    geom_smooth(method = 'lm', formula = y ~ x, se = F)+
    guides(color = F) +
    labs(x = 'Horas estudadas pelo estudante', y = "Desempenho escolar")+
    theme_bw()
)

ggplotly(
  estudante_escola %>% 
    ggplot(aes(x = horas, y = desempenho, color = escola))+
    geom_point()+
    geom_smooth(method = 'lm', formula = y ~ x, se = F)+
    labs(x = 'Horas estudadas pelo estudante', y = "Desempenho escolar")+
    scale_color_viridis_d() +
    theme_bw()
)

base_exemplo <- estudante_escola %>% 
  filter(escola %in% c('1', '2', '3', '4', '5', '6')) %>% 
  mutate(escola = as.numeric(escola))

scatter3d(desempenho ~ horas + escola, 
          data = base_exemplo, fit = 'quadratic')

scatter3d(desempenho ~ horas + escola, groups = factor(base_exemplo$escola),
          data = base_exemplo, fit = 'quadratic', surface = T)

modelo_nulo_hlm2 <- lme(fixed = desempenho ~ 1,
                        random = ~ 1 | escola,
                        data = estudante_escola, 
                        method = 'REML')

summary(modelo_nulo_hlm2)
stderr_nlme(modelo_nulo_hlm2)

modelo_ols_nulo <- lm(formula = desempenho ~ 1, 
                      data = estudante_escola)

summary(modelo_ols_nulo)
lrtest(modelo_ols_nulo, modelo_nulo_hlm2)


data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())


modelo_intercept_hlm2 <- lme(fixed = desempenho ~ horas, 
                             random = ~ 1 | escola, 
                             data = estudante_escola, method = "REML")

summary(modelo_intercept_hlm2)
stderr_nlme(modelo_intercept_hlm2)


data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2), 
           HLM2_Intercept = logLik(modelo_intercept_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2, 
         `HLM2_Intercept` = 3) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45", 'aquamarine4')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

glimpse(estudante_escola)
modelo_intercept_inclin_hlm2 <- lme(fixed = desempenho ~ horas, 
                                    random = ~ horas | escola,
                                    data = estudante_escola,
                                    method = 'REML')

summary(modelo_intercept_inclin_hlm2)
stderr_nlme(modelo_intercept_inclin_hlm2)

data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2), 
           HLM2_Intercept = logLik(modelo_intercept_hlm2),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2, 
         `HLM2_Intercept` = 3,
         `HLM2_Intercept_Inclin_Aleat` = 4) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 4) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45", 'aquamarine4', 'darkorchid')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())


# Modelo final ------------------------------------------------------------

modelo_final_hlm2 <- lme(fixed = desempenho ~ horas + texp + texp:horas,
                         random = ~ horas | escola, 
                         data = estudante_escola, 
                         method = 'REML')

summary(modelo_final_hlm2)
stderr_nlme(modelo_final_hlm2)

data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2), 
           HLM2_Intercept = logLik(modelo_intercept_hlm2),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2),
           HLM2_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2, 
         `HLM2_Intercept` = 3,
         `HLM2_Intercept_Inclin_Aleat` = 4,
         `HLM2_Final` = 5) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 4) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45", 'aquamarine4', 'darkorchid',
                               'darkblue')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

v_final <- data.frame(modelo_final_hlm2$coefficients$random$escola) %>% 
  rename(v00 = 1, v10 = 2)

v_final$escola <- c(1:10)
v_final$escola <- as.factor(v_final$escola)

v_final %>% 
  select(escola, everything()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

random.effects(modelo_final_hlm2) %>% 
  rename(v1j = 2) %>% 
  rownames_to_column("Escola") %>% 
  mutate(color_v1j = ifelse(v1j < 0, yes = "A", no = "B"),
         hjust_v1j = ifelse(v1j > 0, 1.15, -0.15)) %>%
  ggplot(aes(label = format(v1j, digits = 2),
             hjust = hjust_v1j)) +
  geom_bar(aes(x = fct_rev(Escola), y = v1j, fill = color_v1j), 
           stat = 'identity', color = 'black') +
  geom_text(aes(x = Escola, y = 0), size = 4.1, color = "Black") +
  coord_flip() + 
  labs(x = "Escola", y = expression(nu[1][j])) +
  scale_fill_manual(values = c("darkred", "darkgreen"))+
  theme(panel.background = element_rect('white'), 
        panel.border = element_rect(NA), 
        panel.grid = element_line('grey99'),
        legend.position = 'none')


# Incluindo os fitted values do HLM2 no Data frame original --------------------

estudante_escola$hlm2_fitted <- predict(modelo_final_hlm2, estudante_escola)

predict(modelo_final_hlm2, level = 0:1) %>% 
  mutate(escola = gsub("^.*?\\/", "", escola), 
         escola = as.factor(as.numeric(escola)),
         desempenho  = estudante_escola$desempenho,
         etjk = resid(modelo_final_hlm2)) %>% 
  select(escola, desempenho, everything()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                font_size = 12, full_width = T)


# Exemplos de predições ---------------------------------------------------

#Exemplo: Quais os valores previstos de desempenho escolar, para dado
#aluno que estuda na escola "1", sabendo-se que ele estuda 11h semanais,
#e que a escola oferece tempo médio de experiência de seus professores
#igual a 3.6 anos?

predict(object = modelo_final_hlm2, level = 0:1,
        newdata =  data.frame(escola = '1', 
                              horas = 11,
                              texp = 3.6))


# Plotagem dos Valores previstos considerando os efeitos aleatórios de 
# intercepto e inclinação

estudante_escola %>% 
  ggplot(aes(horas, y = hlm2_fitted)) +
  geom_point() +
  geom_smooth(aes(color = escola),method = 'lm', se = F) +
  scale_color_viridis_d()+
  labs(x = "Quantidade de horas de estudo do estudante", 
       y = "Desempenho (Fitted values)", color = "Escola") + 
  theme_bw()


################################################################################
#                       COMPARAÇÃO COM UM MODELO OLS                           #
################################################################################

# Elaborando um modelo OLS (mínimos quadrados ordinários), para efeitos de
# comparação

modelo_ols <- lm(formula = desempenho ~ horas + texp, 
                 data = estudante_escola)

summary(modelo_ols)

# Comparando os LogLiks

data.frame(OLS = logLik(modelo_ols),
           HLM2_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS` = 1,
         `HLM2_Final` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 4) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c('aquamarine4', 'darkorchid')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

# A diferença entre os LL é estatisticamente significante?
lrtest(modelo_final_hlm2, modelo_ols)

# Comparando a aderência dos dois modelos
estudante_escola$ols_fitted <- modelo_ols$fitted.values

estudante_escola %>% 
  ggplot() +
  geom_point(aes(x = desempenho, y = ols_fitted, color = "OLS"), alpha = .3) +
  geom_point(aes(x = desempenho, y = hlm2_fitted, color = "HLM2"), alpha = .3) +
  geom_smooth(aes(x = desempenho, y = desempenho), color = 'darkred', 
              method = 'lm') +
  geom_smooth(aes(x = desempenho, y = ols_fitted), color = 'darkorchid4', 
              method = 'lm', se = F, formula = y ~ splines::bs(x, df=5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y = hlm2_fitted), color = 'aquamarine3', 
              method = 'lm', se = F, formula = y ~ splines::bs(x, df=5),
              size = 1.5)+
  scale_color_manual("Modelos:", values = c('aquamarine3', 'darkorchid4')) +
  labs(x = "Desempenho", y = "Fitted values", color = "Modelos")+
  theme_bw()


################################################################################
#                 COMPARAÇÃO COM UM MODELO OLS COM DUMMIES                     #
################################################################################

# Dummizando (n-1) a variável escola para fazer a função do contexto
estudante_escola_dummies <- dummy_cols(.data = estudante_escola, 
                                       select_columns = 'escola', 
                                       remove_first_dummy = TRUE, 
                                       remove_selected_columns = TRUE)

estudante_escola_dummies %>% 
  select(-hlm2_fitted, -ols_fitted, everything()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', full_width = T, 
                font_size = 12)

# Modelo OLS com dummies
modelo_ols_dummies <- lm(formula = desempenho ~ . - estudante -hlm2_fitted
                         -ols_fitted, data = estudante_escola_dummies)

summary(modelo_ols_dummies)

# Aplicando o procedimento stepwise para retirar as variáveis não
# estatisticamente significante
modelo_ols_dummies_step <- step(object = modelo_ols_dummies, 
                                step = qchisq(p = 0.05, df = 1, 
                                              lower.tail = F))

summary(modelo_ols_dummies_step)

# Comparando os LogLiks

data.frame(OLS = logLik(modelo_ols),
           OLS_Dummies_Step = logLik(modelo_ols_dummies_step),
           HLM2_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS` = 1,
         `OLS_Dummies_Step` = 2,
         `HLM2_Final` = 3) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 4) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c('aquamarine4', 'darkorchid', 'darkgray')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

# A diferença entre os LL é estatisticamente significante?
lrtest(modelo_final_hlm2, modelo_ols_dummies_step)

export_summs(modelo_ols_dummies_step, modelo_final_hlm2, 
             model.names = c("OLS com dummies", "HLM2"))

# Comparando a aderência dos dois modelos
estudante_escola$ols_dummies_step <- modelo_ols_dummies_step$fitted.values

estudante_escola %>% 
  ggplot() +
  geom_point(aes(x = desempenho, y = ols_dummies_step, color = "OLS c/ Dummies"),
             alpha = .3) +
  geom_point(aes(x = desempenho, y = hlm2_fitted, color = "HLM2"), alpha = .3) +
  geom_smooth(aes(x = desempenho, y = desempenho), color = 'darkred', 
              method = 'lm') +
  geom_smooth(aes(x = desempenho, y = ols_dummies_step), color = 'darkorchid4', 
              method = 'lm', se = F, formula = y ~ splines::bs(x, df=5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y = hlm2_fitted), color = 'aquamarine3', 
              method = 'lm', se = F, formula = y ~ splines::bs(x, df=5),
              size = 1.5)+
  scale_color_manual("Modelos:", values = c('aquamarine3', 'darkorchid4')) +
  labs(x = "Desempenho", y = "Fitted values", color = "Modelos")+
  theme_bw()


data.frame(OLS = logLik(modelo_ols),
           HLM2_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS` = 1,
         `HLM2_Final` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 4) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c('aquamarine4', 'darkorchid')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

#Comparação entre os LLs de todos os modelos estimados neste exemplo
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2),
           OLS = logLik(modelo_ols),
           HLM2_Intercept_Aleat = logLik(modelo_intercept_hlm2),
           OLS_Dummies_step = logLik(modelo_ols_dummies_step),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2),
           HLM2_Modelo_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2,
         `OLS` = 3,
         `HLM2 com Interceptos Aleatórios` = 4,
         `OLS com Dummies e Stepwise` = 5,
         `HLM2 com Interceptos e Inclinações Aleatórios` = 6,
         `HLM2 Modelo Final` = 7) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.2, color = "white", size = 5) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_viridis_d() +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())
