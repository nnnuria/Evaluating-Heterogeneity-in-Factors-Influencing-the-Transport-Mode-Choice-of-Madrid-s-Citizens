This project investigates the key determinants of transport mode choice among Madrid cit-
izens using Classification and Regression Trees and weighted average SHAP (SHapley Addit-
ive exPlanations) values while addressing the heterogeneous nature of urban travel behaviour
by applying advanced clustering techniques. The aim is to support urban planners in design-
ing sustainable and inclusive mobility policies. We propose a novel analytical framework
comparing suitability of Latent Class Clustering and an enhanced K-Prototypes algorithm
which incorporates the Eskin measure to handle mixed-type data. These methods are eval-
uated using a dataset consisting of socio-demographic, trip-specific, and weather-related
variables. Model selection criteria such as AIC, silhouette scores, and cluster stability via
bootstrapping are employed to determine the optimal number of clusters. We further assess
the quality of clustering through a classification-based approach using LightGBM. Our find-
ings suggest that Latent Class Clustering outperforms K-Prototypes in extracting meaning-
ful consumer segments. Moreover, across all identified clusters, trip distance, vehicle access,
and possession of a transport card consistently influence mode choice. The SHAP analysis
highlighted the additional impact of weather and rush hours. Thus, the results underscore
the importance of tailoring transport policies to different population groups. Our recom-
mendations advocate for targeted discounts, enhanced safety for children and women, and
better seat availability during peak hours to promote sustainable and inclusive transport in
Madrid.
