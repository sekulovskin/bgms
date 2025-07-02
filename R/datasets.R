#' Post-traumatic stress disorder symptoms of Wenchuan earthquake survivors
#'
#' A data set containing items measuring symptoms of posttraumatic stress
#' disorder (PTSD) \insertCite{McNallyEtAl_2015}{bgms}. Participants were 362
#' Chinese adults who survived the Wenchuan earthquake and lost at least one
#' child in the disaster. PTSD symptoms were reported using the civilian version
#' of the Posttraumatic Checklist, which consists of 17 items, each assessing
#' one of the DSM-IV symptoms of PTSD. Participants rated each item on a
#' five-point scale ranging from  "not at all" to "extremely" to indicate how
#' much the symptom bothered them in the past month.
#'
#' @format A matrix with 362 rows and 17 columns:
#' \describe{
#'   \item{intrusion}{Repeated, disturbing memories, thoughts, or images of a
#'   stressful experience from the past?}
#'   \item{dreams}{Repeated, disturbing dreams of a stressful experience from
#'   the past?}
#'   \item{flash}{Suddenly acting or feeling as if a stressful experience were
#'   happening again (as if you were reliving it)?}
#'   \item{upset}{Feeling very upset when something reminded you of a stressful
#'   experience from the past?}
#'   \item{physior}{Having physical reactions (e.g., heart pounding, trouble
#'   breathing, sweating) when something reminded you of a stressful experience
#'   from the past?}
#'   \item{avoidth}{Avoiding thinking about or talking about a stressful
#'   experience from the past or avoiding having feelings related to it?}
#'   \item{avoidact}{Avoiding activities or situations because they reminded you
#'   of a stressful experience from the past?}
#'   \item{amnesia}{Trouble remembering important parts of a stressful
#'   experience from the past?}
#'   \item{lossint}{Loss of interest in activities that you used to enjoy?}
#'   \item{distant}{Feeling distant or cut off from other people?}
#'   \item{numb}{Feeling emotionally numb or being unable to have loving
#'   feelings for those close to you?}
#'   \item{future}{Feeling as if your future will somehow be cut short?}
#'   \item{sleep}{Trouble falling or staying asleep?}
#'   \item{anger}{Feeling irritable or having angry outbursts?}
#'   \item{concen}{Having difficulty concentrating?}
#'   \item{hyper}{Being "super-alert" or watchful or on guard?}
#'   \item{startle}{Feeling jumpy or easily startled?}
#' }
#'
#' @source \url{http://psychosystems.org/wp-content/uploads/2014/10/Wenchuan.csv}
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @name Wenchuan
#'
#' @usage data("Wenchuan")
#'
#' @references
#'
#' \insertAllCited{}
#'
NULL

#' ADHD symptoms of children aged 6-8 years
#'
#' A data set containing items measuring symptoms of attention deficit
#' hyperactivity disorder (ADHD) \insertCite{Silk_2019_ADHD}{bgms}
#' Participants are a subsample from a community cohort, the Children's Attention
#' Project (CAP); 355 children, aged 6-8 years, of which 146 were diagnosed with
#' ADHD and 209 were not. ADHD symptoms were assessed during structured
#' face-to-face interviews with the children's parents. The NIMH Diagnostic
#' Interview Schedule for Children IV-DISC-iV \insertCite{Shaffer_2000_nimh}{bgms}
#' contains a total of 18 symptoms, of which 9 Inattentive and 9
#' Hyperactive/Impulsive. Each symptom was scored as either present (1) or
#' absent (0), based on the parents' statements.
#'
#' @format A matrix with 355 rows and 19 columns. (I) and (HI) indicate Inattentive
#' and Hyperactive/Impulsive symptoms, respectively.
#' \describe{
#'   \item{group}{1 = ADHD diagnosis, 0 = no ADHD diagnosis}
#'   \item{avoid}{Often avoids, dislikes, or is reluctant to engage in tasks
#'   that require sustained mental effort (I)}
#'   \item{closeatt}{Often fails to give close attention to details or makes
#'   careless mistakes in schoolwork, work, or other activities (I)}
#'   \item{distract}{Is often easily distracted by extraneous stimuli (I)}
#'   \item{forget}{Is often forgetful in daily activities (I)}
#'   \item{instruct}{Often does not follow through on instructions and fails to
#'   finish schoolwork, chores, or duties in the workplace (I)}
#'   \item{listen}{Often does not seem to listen when spoken to directly
#'   (I)}
#'   \item{loses}{Often loses things necessary for tasks or activities (I)}
#'   \item{org}{Often has difficulty organizing tasks and activities (I)}
#'   \item{susatt}{Often has difficulty sustaining attention in tasks or play
#'   activities (I)}
#'   \item{blurts}{Often blurts out answers before questions have been completed
#'   (HI)}
#'   \item{fidget}{Often fidgets with hands or feet or squirms in seat
#'   (HI)}
#'   \item{interrupt}{Often interrupts or intrudes on others (HI)}
#'   \item{motor}{Is often "on the go" or often acts as if "driven by a motor"
#'   (HI)}
#'   \item{quiet}{Often has difficulty playing or engaging in leisure activities
#'   quietly (HI)}
#'   \item{runs}{Often runs about or climbs excessively in situations in which
#'   it is inappropriate (HI)}
#'   \item{seat}{Often leaves seat in classroom or in other situations in which
#'   remaining seated is expected (HI)}
#'   \item{talks}{Often talks excessively (HI)}
#'   \item{turn}{Often has difficulty awaiting turn (HI)}
#' }
#'
#' @source \insertCite{Silk_2019_ADHD;textual}{bgms}. Data retrieved from \url{https://doi.org/10.1371/journal.pone.0211053.s004}.
#' Licensed under the Creative Commons Attribution 4.0 International License (CC-BY 4.0). https://creativecommons.org/licenses/by/4.0/
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @name ADHD
#'
#' @usage data("ADHD")
#'
#' @references
#'
#' \insertAllCited{}
#'
NULL

#' Items measuring boredom proneness
#'
#' A dataset containing items measuring proneness to boredom. The dataset is a
#' subset of the data from \insertCite{Martarelli_2023_Boredom;textual}{bgms},
#' and contains the 8 items from the Short Boredom Proneness Scale (SBPS). Each
#' item was scored on a 7-point Likert scale ranging from 1 ("strongly disagree")
#' to 7 ("strongly agree"). The SBPS was administered either in English (the
#' language of the original scale, by \insertCite{Struk_2015_boredom;textual}{bgms})
#' or French (as translated by \insertCite{Martarelli_2023_Boredom;textual}{bgms}).
#'
#' @format A matrix with 986 rows and 9 columns:
#' \describe{
#'   \item{loose_ends}{I often find myself at “loose ends,” not knowing what to
#'   do.}
#'   \item{entertain}{I find it hard to entertain myself.}
#'   \item{repetitive}{Many things I have to do are repetitive and monotonous.}
#'   \item{stimulation}{It takes more stimulation to get me going than most
#'   people.}
#'   \item{motivated}{I don't feel motivated by most things that I do.}
#'   \item{keep_interest}{In most situations, it is hard for me to find
#'   something to do or see to keep me interested.}
#'   \item{sit_around}{Much of the time, I just sit around doing nothing.}
#'   \item{half_dead_dull}{Unless I am doing something exciting, even dangerous,
#'   I feel half-dead and dull.}
#' }
#'
#' @source \insertCite{Martarelli_2023_Boredom;textual}{bgms}. Data retrieved from \url{https://osf.io/qhux8}.
#' Licensed under the Creative Commons Attribution 4.0 International License (CC-BY 4.0). https://creativecommons.org/licenses/by/4.0/
#'
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @name Boredom
#'
#' @usage data("Boredom")
#'
#' @references
#'
#' \insertAllCited{}
#'
NULL
