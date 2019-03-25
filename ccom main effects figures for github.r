## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----data_prep-----------------------------------------------------------
## setup
library(tidyverse)
library(viridis)
library(broom)
library(broom.mixed)
library(cowplot)
library(patchwork)

ccom <- read_csv("~/Dropbox/PhD/Research/Couples Communication/Food task/Data assembly/Data assembly final/clean/all2.csv")
ccom <- ccom %>%
  mutate(couple = as.factor(couple),
         food = factor(food, levels = c("Radish", "Cookie")),
         pfood = factor(pfood, levels = c("Radish", "Cookie")),
         Total.pinsBin = as.numeric(ccom$Total.pins > 0),
         Total.pins_trunc = ifelse(Total.pins > 0,
                                   Total.pins,
                                   as.numeric(NA)))

couple_ratings <- read_csv("~/Dropbox/PhD/Research/Couples Communication/Food task/Data assembly/Data assembly final/clean/couple_ratings.csv")
couple_ratings$couple <- as.factor(couple_ratings$couple)

couple_data <- couple_ratings %>% select(ID,
                                         couple,
                                         ndep,
                                         couple.debq.emotr,
                                         Mutual.avoidance.,
                                         p.Mutual.avoidance.,
                                         Negative.Reciprocity.,
                                         p.Negative.Reciprocity.,
                                         Positive.Reciprocity.,
                                         p.Positive.Reciprocity.,
                                         Vulnerability.empathy.support.,
                                         p.Vulnerability.empathy.support.) %>%
   mutate_at(vars(-ID, -couple, -ndep),
            funs(scale(-1 / .^(1/2)))) %>%
  gather(variable, value, -ID, -couple, -ndep, -couple.debq.emotr)

## ----figure--------------------------------------------------------------
# data--------------------------------------------------
me_plot_data <- ccom %>%
  select(ID,
         couple,
         ndep,
         food,
         pfood,
         PosAchange,
         PFQ
  ) %>%
  full_join(couple_data %>% 
              filter(variable == "p.Positive.Reciprocity.") %>% 
              spread(variable, value) %>% 
              select(couple, p.Positive.Reciprocity.),
            by = "couple") %>% 
  # reorder food factors
  mutate(food = factor(food, levels = c("Radish", "Cookie")),
         pfood = factor(pfood, levels = c("Radish", "Cookie")),
         ndep = factor(ndep, levels = c("2", "1", "0"))
  ) %>% 
  mutate_at(vars(PosAchange, PFQ),
            funs(
              as.vector(scale(.))
            )
  ) %>%
  gather(variable, value,
         PosAchange,
         PFQ,
         p.Positive.Reciprocity.
  ) %>%
  nest(-variable) %>%
  mutate(se_table = map(.x = data,
                        ~ .x %>%
                          group_by(pfood, food) %>%
                          select(ndep, value) %>%
                          summarise_at(
                            vars(value),
                            funs(
                              ndep = unique(ndep),
                              sd(., na.rm = TRUE),
                              n = sum(!is.na(.)),
                              se=sd/sqrt(n),
                              mean = mean(., na.rm = TRUE)
                                 )
                          )
  )
  )

# theme-------------------------------------------------

food_colors <- viridis(2) %>% 
  set_names(c("Cookie", "Radish"))

theme_set(theme_cowplot())

plot_theme <- list(
  theme(axis.line = element_line(color = 'black'),
        axis.title.y = element_text(
          angle = 0,
          vjust  = 0.5,
          margin = margin(0, 0, 0, 0.5)
        ),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
        plot.subtitle = element_text(face = "bold")
  ),
  scale_y_continuous(position = "right"),
  scale_colour_manual(values = food_colors),
  scale_fill_manual(values = food_colors),
  guides(color = guide_legend(reverse = TRUE)),
  facet_wrap(~ ndep,
             ncol = 1,
             labeller = as_labeller("ndep_labels"),
             strip.position = "bottom",
             scales = "free_y"
  ),
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ),
  coord_flip(),
  # suppress y axis
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        # suppress facet labels
        strip.text.x = element_blank())
)

no_legend <- list(
  theme(legend.position = "none")
)

ndep_labels <- c(
  "0" = "Neither Partner in Cookie Condition",
  "1" = "One Partner in Cookie Condition",
  "2" = "Both Partner in Cookie Condition"
)

# food_colors <- c(
#   "Radish" = "firebrick1",
#   "Cookie" = "chocolate4"
# )

background_colors <- tibble(
  ndep = factor(c(0, 1, 1, 2),
                levels = c("2", "1", "0"))
  ) %>% 
  mutate(
    xmin = c(0.5, 0.5, 1.5,0.5),
    xmax = c(1.5, 1.5, 2.5,1.6),
    fill = c("Radish",
             "Radish",
             "Cookie",
             "Cookie"
    )
  )

background_colors

# cutom legend fill-------------------------------
legend_theme <- theme(
  legend.margin = margin(0, 0, 20, 0),
  legend.justification="center"
)

legend_fill_1 <- ggplot(data = filter(background_colors,
                                      fill == "Cookie"),
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            )
) +
  geom_rect(alpha = 0.2) +
  scale_fill_manual(values = food_colors) +
  labs(fill = "Partners'\nFood")  +
  legend_theme

legend_fill_1 <- get_legend(legend_fill_1)

legend_fill_2 <- ggplot(data = background_colors,
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            )
) +
  geom_rect(alpha = 0.2) +
  scale_fill_manual(values = food_colors) +
  labs(fill = "")  +
  legend_theme

legend_fill_2 <- get_legend(legend_fill_2)

legend_fill_3 <- ggplot(data = filter(background_colors,
                                      fill == "Radish"),
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            )
) +
  geom_rect(alpha = 0.2) +
  scale_fill_manual(values = food_colors) +
  labs(fill = "") +
  legend_theme

legend_fill_3 <- get_legend(legend_fill_3)

legend_fill <- plot_grid(
  legend_fill_1,
  legend_fill_2,
  legend_fill_3,
  ncol = 1, align = "hv"
  )

# cutom legend color-------------------------------
legend_color_1 <-  ggplot(data = me_plot_data %>%
                filter(variable == "p.Positive.Reciprocity.") %>%
                unnest(data) %>% 
                  filter(food == "Cookie"),
              aes(
                x = pfood,
                y = value,
                color = food
              )
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = food_colors) +
  labs(color = "Participants'\nFood") +
  legend_theme

legend_color_1 <- get_legend(legend_color_1)

legend_color_2 <- ggplot(data = me_plot_data %>%
                filter(variable == "p.Positive.Reciprocity.") %>%
                unnest(data),
              aes(
                x = pfood,
                y = value,
                color = food
              )
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = food_colors) +
  labs(color = "") +
  legend_theme

legend_color_2 <- get_legend(legend_color_2)

legend_color_3 <- ggplot(data = me_plot_data %>%
                filter(variable == "p.Positive.Reciprocity.") %>%
                unnest(data) %>% 
                  filter(food == "Radish"),
              aes(
                x = pfood,
                y = value,
                color = food
              )
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = food_colors) +
  labs(color = "") +
  legend_theme

legend_color_3 <- get_legend(legend_color_3)

legend_color <- plot_grid(
  legend_color_1,
  legend_color_2,
  legend_color_3,
  ncol = 1, align = "hv")

# pos rec error bars------------------------------------
# four error bars, one for each food pfood crossing
pos_rec_se_table1 <- me_plot_data %>%
    filter(variable == "p.Positive.Reciprocity.") %>%
    unnest(se_table) %>% 
  add_column(width = 0.5)

# 3 error bars, one for each ndep
pos_rec_se_table2 <- me_plot_data %>%
    filter(variable == "p.Positive.Reciprocity.") %>%
    unnest(se_table) %>%
  mutate(pfood = c(1, 1.5, 1.5, 1)) %>%
  distinct(ndep, .keep_all = TRUE) %>% 
  add_column(width = c(0.8, 1.6, 0.8))

# pos rec plot------------------------------------------
me_plot_posrec <- ggplot(data = pos_rec_se_table2,
                         aes(
                           x = pfood,
                           y = value
                         )
) +
  geom_rect(data = background_colors,
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            ),
            alpha = 0.2
  ) +
  geom_jitter(data = me_plot_data %>%
                filter(variable == "p.Positive.Reciprocity.") %>%
                unnest(data),
              aes(color = food),
              size = 3,
              width = 0.25,
              height = 0.25
  ) +
  labs(x = "Partners' Food",
       y = "Positive Reciprocity in\nAppreciation Discussion (SD)",
       color = "Participants' Food",
       fill = "Partners' Food") +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se,
      x = pfood,
      width = width
      ),
    colour = "black",
    inherit.aes = FALSE
    ) +
  plot_theme +
  geom_point(
    shape=21,
    color = "black",
    fill = "white",
    aes(y = mean),
    size = 3
    )


me_plot_posrec 

# posaff plot------------------------------------------

me_plot_posaff <- ggplot(data = me_plot_data %>%
                           filter(variable == "PosAchange") %>%
                           unnest(se_table),
                         aes(
                           x = pfood,
                           y= value
                         )
) +
  geom_rect(data = background_colors,
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            ),
            alpha = 0.2
  ) +
  geom_jitter(data = me_plot_data %>%
                filter(variable == "PosAchange") %>%
                unnest(data),
              aes(color = food),
              size = 3,
              width = 0.25,
              height = 0.25
  ) +
  labs(
    x = "Partners' Food",
    y = "Change in Positive Affect\nin Disagreement Discussion (SD)"
  ) +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se,
      x = pfood
      ),
    colour = "black",
    width = .8,
    inherit.aes = FALSE
  ) +
  plot_theme +
  geom_point(
    shape=21,
    color = "black",
    fill = "white",
    aes(y = mean),
    size = 3
  )

me_plot_posaff

# pfq plot------------------------------------------

me_plot_pfq <- ggplot(data = me_plot_data %>%
                        filter(variable == "PFQ") %>%
                        unnest(se_table),
                      aes(
                        x = pfood,
                        y= value
                      )
) +
  geom_rect(data = background_colors,
            aes(
              x = NULL,
              y = NULL,
              xmin = xmin, xmax = xmax, 
              ymin = -Inf, ymax = Inf,
              fill = fill
            ),
            alpha = 0.2
  ) +
  geom_jitter(data = me_plot_data %>%
                filter(variable == "PFQ") %>%
                unnest(data),
              aes(color = food),
              size = 3,
              width = 0.25,
              height = 0.25
  ) +
  geom_errorbar(
    aes(
    ymin = mean - se,
    ymax = mean + se,
    x = pfood
  ),
  colour = "black",
  width = .8,
  inherit.aes = FALSE
  ) +
  plot_theme +
  geom_point(
    shape=21,
    color = "black",
    fill = "white",
    aes(y = mean),
    size = 3
    ) +
  labs(y = "Positive Feelings for Partner\nAfter Disagreement Discussion (SD)")

me_plot_pfq

# group plot---------------------------------------------
me_legend <- get_legend(me_plot_posrec +
                          theme(legend.justification="center"))

me_plot_1 <- me_plot_posaff + no_legend + labs(subtitle = "A") +
    me_plot_pfq + no_y_axis + no_legend + labs(subtitle = "B") +
    me_plot_posrec + no_y_axis + no_legend + labs(subtitle = "C") +
  plot_layout(nrow = 1)
  

me_plot_2 <- legend_color + legend_fill + me_plot_1 + plot_layout(nrow = 1, widths = c(1, 1, 10))



save_plot("/Users/Ben/Dropbox/PhD/Research/Couples Communication/Paper draft/main-effects-plot.pdf",
          me_plot_2,
          nrow = 1,
          ncol = 3,
          base_height = 4)

