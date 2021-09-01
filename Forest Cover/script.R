library(tidyverse)


forest <- read.csv('covtype.csv', header = TRUE)

wilderness.count <- forest %>% group_by(Wilderness_Area1, Wilderness_Area2, 
                                        Wilderness_Area3, Wilderness_Area4,
                                        Cover_Type) %>%
  summarise(counts = n())
wilderness.count


ggplot(data = wilderness.count %>% filter(Wilderness_Area1 == 1), 
       aes(x = factor(Cover_Type), y = counts)) +
  geom_bar(stat = 'identity')
