\name{readmission}
\docType{data}
\alias{readmission}

\title{Rehospitalization colorectal cancer}
\usage{data(readmission)}
\description{
 This contains rehospitalization times after surgery in patients diagnosed with colorectal cancer
}
\format{
  This data frame contains the following columns:
  \tabular{ll}{
    id:\tab identificator of each subject (repeated for each recurrence)\cr
    enum:\tab observation number within patient\cr
    t.start:\tab start of interval (0 or previous recurrence time)\cr
    t.stop:\tab recurrence or censoring time\cr
    time:\tab interocurrence or censoring time\cr
    event:\tab rehospitalization status. All event are 1 for each subject excepting last one that it is 0\cr
    chemo:\tab Treatment (chemotherapy) indicator\cr
    sex:\tab Gender\cr
    dukes:\tab Tumour stage under Dukes's classification: A-B, C, or D.\cr
    charlson:\tab Comorbidity Charlson's index, modelled as a time
    dependent covariate and classified into three groups: Index 0, Index 1-2, and Index >= 3\cr
    death:\tab death indicator: Dead = 1; alive = 0\cr
  }
}
\source{
Gonzalez, JR., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M., Cambray, M. and Borras, JM (2005). Sex differences in hospital readmission among colorectal cancer patients. \emph{Journal of Epidemiology and Community Health}, \bold{59}, 6, 506-511.
}

\keyword{datasets}
