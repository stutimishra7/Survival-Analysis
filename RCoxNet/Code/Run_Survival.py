from Run_for_entiredata import *
import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from lifelines.statistics import logrank_test
from matplotlib.backends.backend_pdf import PdfPages


pdf_pages = PdfPages('tcga_OV_Survival_PanCancer Atlas_WITH_RWR_ON_TEST_DATA.pdf')



median_tmb = surv_data['PI'].median()
print(median_tmb)

surv_data['risk_group'] = (surv_data['PI'] >= median_tmb).astype(int)

kmf = KaplanMeierFitter()


high_risk_data = surv_data[surv_data['risk_group'] == 1]
low_risk_data = surv_data[surv_data['risk_group'] == 0]
len(high_risk_data)

kmf.fit(high_risk_data['OS_MONTHS'], event_observed=high_risk_data['OS_EVENT'], label=f'High msi count:{len(high_risk_data)}, msi >= {median_tmb}')
ax = kmf.plot()

kmf.fit(low_risk_data['OS_MONTHS'], event_observed=low_risk_data['OS_EVENT'], label=f'Low msi count:{len(low_risk_data)}, msi <= {median_tmb}')
ax = kmf.plot(ax=ax)

ax.set_title("Survival plot on OV(Test data) with RWR")
ax.set_xlabel("Months")
ax.set_ylabel("Survival Probability")
result = logrank_test(high_risk_data['OS_MONTHS'], low_risk_data['OS_MONTHS'], event_observed_A=high_risk_data['OS_EVENT'], event_observed_B=low_risk_data['OS_EVENT'])
p_value = result.p_value
print(p_value)

plt.text(0.6, 0.7, f'p-value: {p_value:.4f}', transform=ax.transAxes, fontsize=12)


pdf_pages.savefig()
plt.show()

pdf_pages.close()
