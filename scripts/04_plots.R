library(forestplot)
library(dplyr)
library(tibble)

# forest plots of baseline similarity ----

dat_sample1 <- tibble(
  variable = c("Agreeableness", "Consientiousness", "Extraversion", "Neuroticism", "Openness",
               "Attachment Avoidance", "Attachment Anxiety",
               "Responsiveness", "Trust"),
  mean = c(0.191, 0.155, 0.025, 0.069, 0.246,
           0.074, 0.153,
           0.455, 0.236),
  lower = c(0.046, 0.008, -0.123, -0.079, 0.103,
            -0.075, 0.005,
            0.329, 0.089),
  upper = c(0.328, 0.295, 0.171, 0.213, 0.379,
            0.221, 0.295,
            0.566, 0.373),
  similarity = c(".19", ".16", ".03", ".07", ".25",
                 ".07", ".15",
                 ".46", ".24"))

dat_sample1 %>%
  forestplot(labeltext = c(variable, similarity),
             zero = 0)

dat_sample2 <- tibble(
  variable = c("Agreeableness", "Consientiousness", "Extraversion", "Neuroticism", "Openness",
               "Attachment Avoidance", "Attachment Anxiety",
               "Social Support", "Social Provision",
               "Caregiving Proximity", "Caregiving Sensitivity",
               "Caregiving Cooperation", "Compulsive Caregiving",
               "Verbal Aggression", "Collaboration", "Stalemate", "Avoidance-Capitulation"),
  mean = c(0.124, 0.065, -0.001, 0.004, 0.159,
           0.175, 0.195,
           0.084, 0.311,
           0.057, 0.267,
           0.073, 0.043,
           0.421, 0.246, 0.206, 0.219),
  lower = c(-0.028, -0.088, -0.153, -0.148, 0.008,
            0.024, 0.045,
            -0.068, 0.168,
            -0.095, 0.120,
            -0.079, -0.109,
            0.288, 0.098, 0.056, 0.070),
  upper = c(0.271, 0.214, 0.150, 0.155, 0.303,
            0.318, 0.336,
            0.233, 0.442,
            0.207, 0.402,
            0.222, 0.193,
            0.538, 0.383, 0.346, 0.358),
  similarity = c(".12", ".07", "-.00", ".00", ".16",
                 ".18", ".20",
                 ".08", ".31",
                 ".06", ".27",
                 ".07", ".04",
                 ".42", ".25", ".21", ".22"))

dat_sample2 %>%
  forestplot(labeltext = c(variable, similarity),
             zero = 0)

# forest plots of perceived versus actual similarity ----

dat_sample1 <- tibble(
  variable = c("Agreeableness", "", "Consientiousness", "", "Extraversion", "",
               "Neuroticism", "", "Openness", ""),
  mean = c(0.191, 0.330, 0.155, 0.217, 0.025, 0.171, 
           0.069, 0.151, 0.246, 0.375),
  lower = c( 0.046, 0.194, 0.008, 0.074, -0.123, 0.026, 
            -0.079, 0.006, 0.103, 0.243),
  upper = c(0.328, 0.454, 0.295, 0.352, 0.171, 0.309, 
            0.213, 0.290, 0.379, 0.494),
  similarity = c(".19", ".33", ".16", ".28", ".06", ".17", 
                 ".07", ".15", ".25", ".38"))

styles_sample1 <- fpShapesGp(
  lines = list(
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash")
  ),
  box = list(
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c")
  ) 
)

dat_sample1 %>%
  forestplot(labeltext = c(variable, similarity),
             zero = 0, xticks = c(0, 0.5),
             shapes_gp = styles_sample1)

dat_sample2 <- tibble(
  variable = c("Social Support", "", "Social Provision", "",
               "Caregiving Proximity", "", "Caregiving Sensitivity", "",
               "Caregiving Cooperation", "", "Compulsive Caregiving", "",
               "Verbal Aggression", "", "Collaboration", "", 
               "Stalemate", "", "Avoidance-Capitulation", ""),
  mean = c(0.084, 0.475, 0.311, 0.637,
           0.057, 0.365, 0.267, 0.537,
           0.073, 0.346, 0.043, 0.264,
           0.421, 0.597, 0.246, 0.649, 
           0.206, 0.479, 0.219, 0.479),
  lower = c(-0.068, 0.349,  0.168, 0.537,
            -0.095, 0.226,  0.120, 0.420,
            -0.079, 0.205, -0.109, 0.117,
             0.288, 0.490,  0.098, 0.551, 
             0.056, 0.353,  0.070, 0.353),
  upper = c(0.233, 0.585, 0.442, 0.719,
            0.207, 0.489, 0.402, 0.637,
            0.222, 0.473, 0.193, 0.399,
            0.538, 0.686, 0.383, 0.729, 
            0.346, 0.588, 0.358, 0.588),
  similarity = c(".08", ".48", ".31", ".64",
                 ".06", ".37", ".27", ".54",
                 ".07", ".35", ".04", ".26",
                 ".42", ".60", ".25", ".65", 
                 ".21", ".48", ".22", ".48"))

styles_sample2 <- fpShapesGp(
  lines = list(
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash"),
    gpar(col = "#E69F00",
         lty = "solid"),
    gpar(col = "#710c0c",
         lty = "longdash")
  ),
  box = list(
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c"),
    gpar(fill = "#E69F00"),
    gpar(fill = "#710c0c")
  ) 
)

dat_sample2 %>%
  forestplot(labeltext = c(variable, similarity),
             zero = 0,
             shapes_gp = styles_sample2)