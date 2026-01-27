################################################################################
# Trigeminal Sensitivity Analysis - Bormann Study
# Author: Jorn Lotsch
# Description: Analysis of trigeminal sensitivity study data, including
#              preprocessing, translation, categorization, descriptive stats,
#              and tabulations per variable category.
################################################################################

# ======================== #
# Functions
# ======================== #

#' Custom publication-quality ggplot2 theme
#'
#' @return ggplot2 theme object
theme_plot <- function() {
  theme_minimal(base_family = "Libre Franklin") +
    theme(
      plot.title = element_text(
        face = "plain", size = 12, color = "#222222",
        hjust = 0, margin = margin(b = 10)
      ),
      axis.title = element_text(face = "plain", size = 10, color = "#444444"),
      axis.text = element_text(face = "plain", size = 10, color = "#444444"),
      plot.caption = element_text(
        size = 8, color = "#888888",
        hjust = 0, margin = margin(t = 10)
      ),
      panel.grid.major.y = element_line(
        color = "#dddddd", linetype = "dashed", size = 0.3
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks = element_line(color = "#bbbbbb", size = 0.5),
      axis.ticks.length = unit(5, "pt"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = margin(20, 20, 20, 20),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "#222222")
    )
}

#' Scale vector to specified range
#'
#' @param x Numeric vector to scale
#' @param minX Desired minimum value
#' @param maxX Desired maximum value
#' @return Scaled vector in range [minX, maxX]
scaleRange_01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Min-max normalization to [0,1] with fixed boundaries
#'
#' @param x Numeric vector to normalize
#' @param minX Fixed minimum boundary
#' @param maxX Fixed maximum boundary
#' @return Normalized vector
scale01minmax <- function(x, minX, maxX) {
  (x - minX) / (maxX - minX)
}

# ======================== #
# Global variables
# ======================== #

# Define category label mapping
category_labels <- c(
  "1" = "Demographics",
  "2" = "Disorders or health complaints",
  "3" = "COVID-19 history",
  "4" = "Smoking and alcohol use",
  "5" = "Subjective nasal trigeminal chemosensory perception",
  "6" = "Rated olfactory function",
  "7" = "Ratings of nasal irritation and airflow",
  "8" = "Psychophysical measurements"
)

# Variable name translations
translation_map <- c(
# Category 1 - Demographics
  "Alter" = "Age",
  "Geschlecht" = "Gender",
  "Gewicht" = "Weight",
  "Körpergröße" = "Height",

# Category 2 - Disorders or health complaints
  "Allergische Probleme" = "Allergic problems",
  "Infektion der oberen Atemwege" = "Upper respiratory infection",
  "Chronische Sinusitis" = "Chronic sinusitis",
  "OP im HNO-Bereich" = "Surgery in ENT region",
  "chronische Erkrankungen" = "Chronic disease",
  "neurologische Erkrankung" = "Neurological disorder",
  "Ärztliche Vorstellung wegen Problematik der Nasenatmung" = "Medical consultation for nasal breathing problems",
  "Wenn ja erfolgte eine Therapie? - Welche?" = "If yes, was therapy performed (and which one)",
  "Gesichtsschmerzen" = "Facial pain",
  "Wie oft haben Sie Gesichtsschmerzen?" = "How often do you have facial pain",
  "Wie sind die Gesichtsschmerzen?" = "What is the nature of facial pain",
  "Wie oft haben Sie im Monat Migräne?" = "How often do you have migraine per month",
  "Hat sich die Migräne in den letzten 10 Jahren verändert?" = "Has migraine changed over the last 10 years",

# Category 3 - COVID-19 history
  "Waren Sie bereits an Covid erkrankt?" = "Have you had COVID-19?",
  "R1 in %" = "Smell ability before COVID-19",
  "Besteht oder bestand eine Riechminderung nach Covid-19?" = "Is or was there smell reduction after COVID-19?",
  "R2 in %" = "Smell ability immediately after COVID-19",
  "Wie oft Covid-19?" = "How many times have you had COVID-19?",
  "Zeitraum 1" = "Period 1",
  "Riechminderung?" = "Smell reduction 1",
  "Zeitraum 2" = "Period 2",
  "Riechminderung?2" = "Smell reduction 2",
  "Zeitraum 3" = "Period 3",
  "Riechminderung?3" = "Smell reduction 3",

# Category 4 - Smoking and alcohol use
  "Rauchen Sie?" = "Do you smoke?",
  "Wenn Ja: Wie viele Zigaretten am Tag?" = "If yes: how many cigarettes per day?",
  "Wenn ja: Seit wann?" = "If yes: since when?",
  "Waren Sie je Raucher?" = "Have you ever smoked?",
  "Wenn ja: In welchem Zeitraum?" = "If yes: in what time period?",
  "Wenn ja: Wie viele Zigaretten am Tag?5" = "If yes: how many cigarettes then per day?",
  "Trinken Sie Alkohol?" = "Do you drink alcohol?",

# Category 5 - Subjective nasal chemosensory perception (TriFunQ and related)
  "Stechende oder brennende Gerüche wie Rauch, Essig oder Nagellackentferner rufen starke Emotionen in mir hervor" =
    "Pungent or burning odors like smoke, vinegar, or nail polish remover elicit strong emotions in me",
  "Wenn ich einen frischen Minzkaugummi gekaut habe, habe ich das Gefühl besser Luft durch die Nase zu bekommen" =
    "After chewing a fresh mint gum, I feel I can breathe better through my nose",
  "Ich meide kohlensäurehaltige Getränke, weill sie in der Nase z.B. beim Aufstoßen zu sehr brennen" =
    "I avoid carbonated beverages because they burn my nose (e.g., on burping)",
  "Ich schätze meine Nasenatmung als sehr gut ein" =
    "I consider my nasal breathing to be very good",
  "Ich gehe nicht gerne in eine Sauna, weil ich heiße Luft in der Nase als brennend wahrnehme" =
    "I dislike going to saunas because I perceive hot air in my nose as burning",
  "Beim Zwiebelschneiden tränen mir die Augen stark" =
    "My eyes tear strongly when cutting onions",
  "Stechende oder brennende Gerüche lösen bei mir Husten oder Niesen aus" =
    "Pungent or burning odors cause me to cough or sneeze",
  "Ich meide brennende, stechende Gerüche (z.B. Ammoniak o. Chlor)" =
    "I avoid burning or pungent smells (e.g., ammonia or chlorine)",
  "Wenn ich etwas Beißendes oder Stechendes rieche, denke ich voller Panik an Situationen, in denen mir so etwas ähnliches passiert ist" =
    "When I smell something biting or pungent, I panic, remembering similar situations",
  "Im Winter ist mir die kalte Luft in der Nase ausgesprochen unangenehm" =
    "In winter, I find cold air in my nose extremely uncomfortable",
  "Wenn ich Meerrettich esse, finde ich das Beißen in der Nase besonders störend" =
    "When I eat horseradish, I find the burning in my nose especially bothersome",
  "Brennende oder stechende Gerüche können bei mir zu unangenehmen Empfindungen/ Schmerzen im Gesicht führen" =
    "Burning or pungent odors can cause unpleasant sensations or pain in my face",
  "Ich nehme die Kohlensäure in Getränken bewusst und intensiv wahr" =
    "I consciously and intensely perceive carbonation in drinks",
  "Wenn es um leicht kitzelnde oder stechende Gerüche geht, ist meine Nase viel empfindlicher für das Stechende und Beißende, als die Nase anderer Leute" =
    "When it comes to slightly tingling or pungent odors, my nose is much more sensitive than others'",
  "Ich benutze nur Zahnpasta mit sehr mildem Minzgeruch" =
    "I only use toothpaste with a very mild mint scent",
  "Wie oft schneiden Sie im Monat frische Zwiebeln?" =
    "How often do you cut fresh onions per month?",
  "Wenn Sie über das letzte halbe Jahr nachdenken, wie stark tränen Ihnen die Augen beim Zwiebelnschneiden?" =
    "Thinking back over the past six months, how much did your eyes tear when cutting onions?",
  "Hat sich das Tränen Ihrer Augen in den letzten 10 Jahren verändert?" =
    "Has your eye watering while cutting onions changed in the last 10 years?",

# Category 6 - Rated olfactory function
  "R3 in %" = "Current smell ability",
  "Riech- und Schmeckvermögen vermindert" = "Reduced smell and taste ability",
  "Wenn ja, wie begann das Problem" = "If yes, how did the problem start",
  "Wie hat sich das Problem verändert?" = "How has the problem changed",

# Category 7 - Ratings of nasal irritation and airflow
  "Wie empfindlich ist Ihre Nase für stechendes/brennendes?" = "Sensitivity of the nose to stinging/burning stimuli",
  "R23" = "Nasal airflow (both nostrils)",
  "R24" = "Nasal airflow (right nostril)",
  "R25" = "Nasal airflow (left nostril)",

# Category 8 - Psychophysical measurements
  "Lateralisierung (x/20)" = "Lateralization (x/20)",
  "R28" = "AmmoLa intensity",
  "Identifikationstest (x/3)" = "Odor identification (x/3)",
  "CO2-Schwelle" = "CO2 threshold"
)


variables_by_categories <- list(
  Demographics = c("Age", "Gender", "Weight", "Height"),
  Disorders_or_health_complaints = c(
    "Allergic problems",
    "Upper respiratory infection",
    "Chronic sinusitis",
    "Surgery in ENT region",
    "Chronic disease",
    "Neurological disorder",
    "Medical consultation for nasal breathing problems",
    "If yes, was therapy performed (and which one)",
    "Facial pain",
    "How often do you have facial pain",
    "What is the nature of facial pain",
    "How often do you have migraine per month",
    "Has migraine changed over the last 10 years"
  ),
  COVID19_history = c(
    "Have you had COVID-19?",
    "Smell ability before COVID-19",
    "Is or was there smell reduction after COVID-19?",
    "Smell ability immediately after COVID-19",
    "How many times have you had COVID-19?",
    "Period 1",
    "Smell reduction 1",
    "Period 2",
    "Smell reduction 2",
    "Period 3",
    "Smell reduction 3"
  ),
  Smoking_and_alcohol_use = c(
    "Do you smoke?",
    "If yes: how many cigarettes per day?",
    "If yes: since when?",
    "Have you ever smoked?",
    "If yes: in what time period?",
    "If yes: how many cigarettes then per day?",
    "Do you drink alcohol?"
  ),
  Nasal_chemosensory_perception = c(
    "Pungent or burning odors like smoke, vinegar, or nail polish remover elicit strong emotions in me",
    "After chewing a fresh mint gum, I feel I can breathe better through my nose",
    "I avoid carbonated beverages because they burn my nose (e.g., on burping)",
    "I consider my nasal breathing to be very good",
    "I dislike going to saunas because I perceive hot air in my nose as burning",
    "My eyes tear strongly when cutting onions",
    "Pungent or burning odors cause me to cough or sneeze",
    "I avoid burning or pungent smells (e.g., ammonia or chlorine)",
    "When I smell something biting or pungent, I panic, remembering similar situations",
    "In winter, I find cold air in my nose extremely uncomfortable",
    "When I eat horseradish, I find the burning in my nose especially bothersome",
    "Burning or pungent odors can cause unpleasant sensations or pain in my face",
    "I consciously and intensely perceive carbonation in drinks",
    "When it comes to slightly tingling or pungent odors, my nose is much more sensitive than others'",
    "I only use toothpaste with a very mild mint scent",
    "How often do you cut fresh onions per month?",
    "Thinking back over the past six months, how much did your eyes tear when cutting onions?",
    "Has your eye watering while cutting onions changed in the last 10 years?"
  ),
  Rated_olfactory_function = c(
    "Current smell ability",
    "Reduced smell and taste ability",
    "If yes, how did the problem start",
    "How has the problem changed"
  ),
  Nasal_irritation_and_airflow = c(
    "Sensitivity of the nose to stinging/burning stimuli",
    "Nasal airflow (both nostrils)",
    "Nasal airflow (right nostril)",
    "Nasal airflow (left nostril)"
  ),
  Psychophysical_measurements = c(
    "Lateralization (x/20)",
    "AmmoLa intensity",
    "Odor identification (x/3)",
    "CO2 threshold"
  )
)


# Translation list ENT surgeries and nasal breathing problems
condensed_ENT_dictionary <- c(

# Adenoid procedures
  "Adenotomie" = "Adenotomy",
  "Adenotomie Kindesalter" = "Adenotomy",
  "Kindesalter Adenotomie" = "Adenotomy",
  "Adentomie" = "Adenotomy",

# Tonsil procedures
  "Mandeln" = "Tonsillectomy",
  "Mandeln OP" = "Tonsillectomy",
  "Mandel OP" = "Tonsillectomy",
  "Mandel Op" = "Tonsillectomy",
  "Mandel-OP" = "Tonsillectomy",
  "Kindesalter Mandeln" = "Tonsillectomy",
  "Tonsillektomie" = "Tonsillectomy",

# Partial tonsil removal
  "Tonsillotomie" = "Tonsillotomy",
  "Tonsillotomie Kindheit" = "Tonsillotomy",
  "Kindesalter Tonsillotomie" = "Tonsillotomy",
  "Tonsillotomie Kindesalter" = "Tonsillotomy",
  "Tonsillotomie jugend" = "Tonsillotomy",
  "er Tonsillotomie" = "Tonsillotomy",

# Polyp surgery (nasal)
  "Polypen" = "Nasal polyps surgery",
  "Kleinkind Polypen" = "Adenoidectomy",
  "Kindheit Polypen" = "Adenoidectomy",
  "Kindesalter Polypen" = "Adenoidectomy",
  "CRS Polypenentfernung 2x" = "Nasal polyps surgery",
  "Polypen-Op" = "Nasal polyps surgery",
  "OP Polypen Stirnhöhle" = "Sinus surgery",

# Nasal/sinus surgeries and trauma
  "Nasenfraktur" = "Nasal fracture surgery",
  "Nasen # OP" = "Nasal fracture surgery",
  "Wucherungen Nase" = "Nasal polyps surgery",
  "NNH" = "Sinus surgery",
  "NNH Op" = "Sinus surgery",
  "NNH links" = "Sinus surgery",

# Nasal septum procedures
  "NSW" = "Nasal septoplasty",
  "NSW OP" = "Nasal septoplasty",
  "OP NSW" = "Nasal septoplasty",
  "OP Nasenscheidewand" = "Nasal septoplasty",
  "Nasenscheidewand" = "Nasal septoplasty",
  "Nasenscheidewand OP" = "Nasal septoplasty",
  "Begradigung Nasenscheidewand" = "Nasal septoplasty",
  "Begradigung NSW" = "Nasal septoplasty",
  "Septum OP" = "Nasal septoplasty",
  "Septumkorrektur" = "Nasal septoplasty",
  "Septumplastik" = "Nasal septoplasty",
  "Septumplastik li." = "Nasal septoplasty",
  "Septumdeviation" = "Nasal septoplasty",

# Turbinate procedures
  "Nasenmuscheln Verkleinerung" = "turbinate surgery",
  "Reduktion Nasenmuscheln" = "turbinate surgery",
  "Verkleinerung Nasenmuscheln" = "turbinate surgery",
  "Nasenmuschelverkleinerung" = "turbinate surgery",
  "Laserung Nasenmuscheln" = "turbinate surgery",
  "Lasertherapie" = "Nasal mucosa surgery",
  "Laserung Nasenmuscheln" = "turbinate surgery",
  "Verödung Nasenmuschel" = "turbinate surgery",
  "Verödung Nasenmuscheln" = "turbinate surgery",
  "Verödung" = "Nasal mucosa surgery",
  "Stromwellentherapie" = "Radio wave therapy",
  "Stromwellentherapie zum Abschwellen der Nasenmuscheln" = "Radio wave therapy",

# Nose correction / reconstruction
  "Begradigung Nase" = "Nasal septoplasty",

# Ear surgeries
  "Ohr anlegen" = "Otoplasty",
  "Otoplastik" = "Otoplasty",
  "Ohr li. Implantat" = "Ear implant",
  "Knorpelplastik Ohr" = "Otoplasty",
  "Austausch Hämmerchen Ohr" = "Malleus replacement (middle ear)",

# Tympanic membrane surgeries
  "Trommelfell" = "Tympanic membrane surgery",
  "Trommelfell 2x" = "Tympanic membrane surgery",
  "Trommelfell-OP" = "Tympanic membrane surgery",
  "Trommelfellimplantat" = "Tympanic membrane surgery",
  "Trommelfellverschluss" = "Tympanic membrane surgery",
  "Trommelfellschnitt" = "Tympanic membrane surgery",

# Ear trauma
  "Knalltrauma" = "Acoustic trauma",

# Vocal cord / larynx
  "Stimmband OP" = "Vocal cord surgery",
  "Stimmknötchen" = "Vocal cord surgery",

# Oncological surgeries
  "Mundbodenkarzinom" = "Oral carcinoma surgery",
  "Parotisentfernung" = "Parotidectomy",
  "Karzinom Hals" = "Neck carcinoma surgery",
  "Karzinom Rachen li." = "pharyngeal carcinoma surgery",
  "Medulläres Schildrüsenkarzinom" = "Thyroid carcinoma surgery",
  "OP Halsbereich" = "Neck surgery",
  "Abszess Rachen" = "Pharyngeal abscess surgery",
  "Trauma Hals" = "Neck trauma surgery",

# Sinus and jaw procedures
  "Fensterung" = "Fenestration",
  "Fensterung Kieferhöhle" = "Caldwell-Luc procedure",
  "Fensterung Kieferhöhle li." = "Caldwell-Luc procedure",
  "Kiefernzyste" = "Jaw cyst removal",

# Thyroid procedures
  "Schilddrüse" = "Thyroid surgery",
  "Thyreodektomie" = "Thyroid surgery",

# General surgeries
  "OP" = "Operation",
  "Op geplant" = "Operation planned",
  "Trauma Hals" = "Neck trauma surgery",

# Medications and therapies
  "Antibiotika" = "Antibiotics",
  "AB" = "Antibiotics",
  "Cortisonspray" = "Cortisone spray",
  "Cortison-NS" = "Cortisone spray",
  "Cortison-Nasenspray" = "Cortisone spray",
  "Kortison-NS" = "Cortisone spray",
  "Cortison NS" = "Cortisone spray",
  "Cortison NS, Desensibilisierung" = "Cortisone nasal spray and desensitization",
  "Allergie Nasenspray Mometa" = "Cortisone spray",
  "Mometason" = "Cortisone spray",

# Immunotherapy or allergy treatments
  "Desensibilisierung" = "Hyposensitization",
  "Hyposensibilisierung" = "Hyposensitization",
  "Hyposensibilisierung, NT, Tabletten" = "Hyposensitization, nasal therapy, tablets",
  "Antiallergika" = "Antiallergics",
  "Cetirizin" = "Antiallergics",

# Respiratory therapies
  "Atemgerät Schlaf" = "Obstructive Sleep Apnea Syndrome",
  "Asthmadiagnostik" = "Asthma diagnostics",
  "Bedarfsspray Asthma" = "Asthma rescue spray",
  "Corticoidinhalation" = "Corticoid inhalation",

# Nasal remedies
  "Nasenspray" = "Cortisone spray",
  "Nasensalbe" = "Ointment",
  "Meersalzinhalation" = "Sea salt inhalation",

# Miscellaneous
  "Regelmäßige Untersuchungen" = "Regular checkups",
  "Begradigung NSW" = "Nasal septal surgery",
  "Salbe" = "Ointment",
  "Creme" = "Ointment",
  "Doprident" = "Doprident (medication)",
  "Fluimveil" = "acetylcysteine",
  "Rotlichthterapie" = "Red light therapy",
  "Akupunktur" = "Acupuncture",
  "Stirnhöhlenvereiderung?" = "Sinusitis",
  "Pansinusitis" = "Sinusitis",
  "Hals" = "Neck surgery",
  "Verödung Nasenscheidewand" = "Nasal septal surgery",
  "Erweiterung Durchgänge NNH" = "Sinus surgery",
  "Verödung Nase" = "Nasal mucosa surgery",
  "Paukenröhrchen" = "Tympanostomy tube insertion",
  "Paukenröhrchen OP" = "Tympanostomy tube insertion",
  "Knochewnimplantat (Mittelohrschwerhörigkeit)" = "Middle ear surgery",
  "NSW Fensterung" = "Nasal septal surgery",
  "Siebbeinhöhle" = "Sinus surgery",
  "Abtragung Nasenschleimhaut" = "Nasal mucosa surgery",
  "Septumperforation" = "Nasal septal surgery",

# Special cases
  "j" = "Unspecified",
  "n" = "No therapy"
)




