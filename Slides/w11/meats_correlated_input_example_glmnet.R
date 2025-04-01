### meats data example with glmnet
### demonstrates visualizing highly correlated inputs
### demonstrates applying preprocessing to deal with asymmetric distributions
### demonstrates using glmnet (elastic net) to check if the
### input correlation causes problems!!

library(tidyverse)

data('meats', package = 'modeldata')

### the meats data are available within memory!!!
ls()

meats %>% glimpse()

### lots of inputs!!!! several outputs

### look at the histograms of the inputs!

meats %>% 
  select(starts_with('x')) %>% 
  dim()

meats %>% 
  select(starts_with('x')) %>% 
  pivot_longer(everything()) %>% 
  ggplot(mapping = aes(x = value)) +
  geom_histogram(bins = 11) +
  facet_wrap(~name, scales = 'free') +
  theme_bw()

### boxplot to check ranges/scales of the inputs
meats %>% 
  select(starts_with('x')) %>% 
  pivot_longer(everything()) %>% 
  ggplot(mapping = aes(x = name, y = value)) +
  geom_boxplot(fill = 'grey') +
  theme_bw()

### use horizontal boxplots to help with long names
meats %>% 
  select(starts_with('x')) %>% 
  pivot_longer(everything()) %>% 
  ggplot(mapping = aes(y = name, x = value)) +
  geom_boxplot(fill = 'grey') +
  theme_bw()

### or use some string manipulation to pull out the digits!
### the stringr package is loaded with tidyverse!!!
### need to set the group aesthetic when boxplots are grouped
### by numeric variable values!!!!
meats %>% 
  select(starts_with('x')) %>% 
  pivot_longer(everything()) %>% 
  mutate(input_id = stringr::str_extract(name, '\\d+')) %>% 
  mutate(input_id = as.numeric(input_id)) %>% 
  ggplot(mapping = aes(x = input_id, y = value)) +
  geom_boxplot(fill = 'grey',
               mapping = aes(group = input_id)) +
  theme_bw()

### the inputs look slightly skewed, what if we log-transform them?
### can only do this if the inputs are positive and non-zero
### can ADD a number to make all inputs positive if needed
meats %>% 
  select(starts_with('x')) %>% 
  pivot_longer(everything()) %>% 
  ggplot(mapping = aes(x = log(value))) +
  geom_histogram(bins = 11) +
  facet_wrap(~name, scales = 'free') +
  theme_bw()

### "helps" remove skew for some others are still asymmetric
### but nothing extreme

### need to check the correlation structure of the inputs
### use corrplot!!!
meats %>% 
  select(starts_with('x')) %>% 
  cor() %>% 
  corrplot::corrplot()

### use upper triangular to simplify
meats %>% 
  select(starts_with('x')) %>% 
  cor() %>% 
  corrplot::corrplot(type = 'upper')

### I like squares instead of circles
meats %>% 
  select(starts_with('x')) %>% 
  cor() %>% 
  corrplot::corrplot(type = 'upper', method = 'square')

### or use color to easily remove the "grid"
### the squares are all the same size now
meats %>% 
  select(starts_with('x')) %>% 
  cor() %>% 
  corrplot::corrplot(type = 'upper', method = 'color')

### ALOT of input correlation!!!

### there are several outputs, but let's just focus on one of them
meats %>% 
  select(-starts_with('x'))

### we'll focus on protein
meats %>% 
  ggplot(mapping = aes(x = protein)) +
  geom_histogram(bins = 11) +
  theme_bw()

### the output MARGINAL distribution looks odd but it's NOT highly skewed
### however, I bet protein can't be negative...so let's log transform it!
meats %>% 
  ggplot(mapping = aes(x = log(protein))) +
  geom_histogram(bins = 11) +
  theme_bw()

### the output marginal is skewed but...that's ok!!! 
### what matters is the distribution of the RESIDUALS in linear models!!!!
### not the distribution of the MARGINAL output!!!

### need to see the relationship between the response and all inputs
meats %>% 
  pivot_longer(starts_with('x')) %>% 
  ggplot(mapping = aes(x = value, y = protein)) +
  geom_point(alpha = 0.33) +
  facet_wrap(~name) +
  theme_bw()

### look at the relationships of the log inputs to log protein
meats %>% 
  pivot_longer(starts_with('x')) %>% 
  ggplot(mapping = aes(x = log(value), y = log(protein))) +
  geom_point(alpha = 0.33) +
  facet_wrap(~name) +
  theme_bw()

### we have seen that the inputs are slightly skewed
### all inputs are positive so we can easily apply log transformation
### all inputs are VERY correlated

### the output is not highly skewed, but we should consider log transforming
### the output so we can never predict a negative number!

### let's immediately start with regularized non-Bayesian linear models
### in this example (though we should fit the non-Bayesian MLE first to
### help explore!!!)

### should we use ridge or lasso???
### not sure? that's ok!!! let ELASTIC NET figure it out for us!!!!

### create the training set which consists of only inputs and the LOG
### transformed protein
df <- meats %>% 
  mutate(y = log(protein)) %>% 
  select(starts_with('x'), y)

df %>% names()

### use caret to manage everything

library(caret)

### use 5-fold for simplicity
my_ctrl <- trainControl(method = 'cv', number = 5)

my_metric <- 'RMSE'

### use linear additive features for ALL inputs!!!
### we should ALWAYS standardize inputs when using regularized models!!!

set.seed(81231)
fit_enet_add_default <- train( y ~ ., data = df,
                               method = 'glmnet',
                               metric = my_metric,
                               preProcess = c("center", "scale"),
                               trControl = my_ctrl)

fit_enet_add_default

fit_enet_add_default$bestTune

plot(varImp(fit_enet_add_default))

### try a custom grid
my_lambda_grid <- exp(seq(log(min(fit_enet_add_default$results$lambda)),
                          log(max(fit_enet_add_default$results$lambda)),
                          length.out = 15))

enet_grid <- expand.grid(alpha = seq(0.1, 0.9, by = 0.1),
                         lambda = my_lambda_grid)

enet_grid %>% nrow()

set.seed(81231)
fit_enet_add_tune <- train( y ~ ., data = df,
                            method = 'glmnet',
                            metric = my_metric,
                            tuneGrid = enet_grid,
                            preProcess = c("center", "scale"),
                            trControl = my_ctrl)

plot( fit_enet_add_tune, xTrans = log )

fit_enet_add_tune$bestTune

### should we consider a log-transformed input???
### let's try a more general BoxCox!!!
set.seed(81231)
fit_enet_bc_add_default <- train( y ~ ., data = df,
                                  method = 'glmnet',
                                  metric = my_metric,
                                  preProcess = c("BoxCox", "center", "scale"),
                                  trControl = my_ctrl)

fit_enet_bc_add_default
