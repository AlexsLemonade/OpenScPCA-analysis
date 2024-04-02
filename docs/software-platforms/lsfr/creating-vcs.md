# Creating virtual computers with Lightsail for Research

While working on an [analysis module](../../contributing-to-analyses/analysis-modules/index.md), you can develop locally or on a virtual computer with Lightsail for Research (LSfR).

Virtual computers provide access to a set amount of virtual CPUs, memory, and storage through a web browser rather than a physical machine.
Using a virtual computer allows you to develop and run analyses that require more CPUs and memory than is available on your local machine.

_Note that you will need access to AWS to use LSfR._
See [getting access to AWS](../../getting-started/accessing-resources/index.md#getting-access-to-aws) to learn how to get set up with AWS.

- The virtual computers provided through LSfR use an Ubuntu operating system.
- During setup, you will choose your desired configuration, including the amount of virtual CPUs, memory, and storage space.
    - You also have the option to add [additional storage to your virtual computer](STUB-LINK for volumes).
- All virtual computers will have a set of pre-installed applications needed for working with OpenScPCA, including:
    - Git and GitKraken
    - The AWS command line interface
    - R
    - Conda

!!! note "More on virtual computers"

    Read more about [virtual computers with Lightsail for Research](https://docs.aws.amazon.com/lightsail-for-research/latest/ug/computers.html).

## How to create a virtual computer

Follow the below steps to create a virtual computer to use with LSfR.

1. Navigate to the [access portal URL from when you set up your user in IAM Identity center](../aws/index.md#joining-iam-identity-center).
This link should bring you to the AWS Console and home page.

    Once here, be sure that you are in the `us-east-2` region, by selecting the drop-down menu in the tool bar.

    <figure markdown="span">
        ![Select region](../../img/creating-vcs-1.png){width="400"}
    </figure>

1. You will need to use the AWS Service Catalog to create virtual computers.
Open the AWS Service Catalog by using the search bar and typing, `service catalog`.

    <figure markdown="span">
        ![Select service catalog](../../img/creating-vcs-2.png){width="600"}
    </figure>

1. Select `LightsailInstance` in the product list, and click `Launch product`.

    <figure markdown="span">
        ![Launch product](../../img/creating-vcs-3.png){width="600"}
    </figure>

1. You will then choose the name and configurations for your virtual computer.

    - Start by providing a `Provisioned product name`.
    <!--TODO Do we want to provide guidance on names?-->
    - Pick an `Application` to use from the drop-down menu.
    This installs any additional applications along with the pre-installed applications.
        - If you are planning to develop in R, we recommend choosing `Rstudio`, which includes both RStudio and VSCodium.
        - If you are planning to develop in Python, we recommend choosing `VSCodium`, which includes VSCodium, but does not include Rstudio.
        - For information on all application options, see [the LSfR documentation on Applications](https://docs.aws.amazon.com/lightsail-for-research/latest/ug/blueprints-plans.html).
    - Choose any of the availability zones.
    - Name your instance, be sure to make sure this is different than the name chose for your product.
    - Pick the size for your instance.
        - Use this [table describing the total vCPUs, memory, and storage space](https://docs.aws.amazon.com/lightsail-for-research/latest/ug/blueprints-plans.html#plans) to choose the most appropriate instance plan.
        - Before you request a plan with GPUs, you will need to put in a [request with the Data Lab team](../../getting-started/accessing-resources/getting-access-to-compute.md#gpu-instance-access).
    - Set the `ShutdownIdlePercent` to 5 and the `ShutdownTimePeriod` to 10.
    This means if the instance is using less than 5% of the total available CPUs or inactive for more than 10 minutes, the instance will temporarily shut down.
    This _does not_ terminate your instance, it temporarily stops any idle instances until you are ready to resume your work.

    <figure markdown="span">
        ![Configure instance](../../img/creating-vcs-4.png){width="600"}
    </figure>

1. Once you have configured your instance, click `Launch product`.
You have now created a virtual computer!

## How to access a virtual computer

Once you have created your virtual computer, you will be able to access and launch the computer on a web browser.
Note that it can take 5-10 minutes to create the instance.

1. To access a virtual computer, you will need to navigate to Lightsail for Research.
Search for Lightsail for Research using the search bar and then select the product.

    <figure markdown="span">
        ![Search for Lightsail](../../img/creating-vcs-5.png){width="600"}
    </figure>

1. This will take you directly to a page that lists your virtual computers.
You should see the virtual computer you created with the instance name that you provided during set up.

    <figure markdown="span">
        ![List of computers](../../img/creating-vcs-6.png){width="600"}
    </figure>

1. To launch the computer, click `Start computer`.
When the computer is ready, use the drop-down menu in the lower right-hand side to select `Access operating system`.

    <figure markdown="span">
        ![Start computer](../../img/creating-vcs-7.png){width="600"}
    </figure>

    You can also use `Launch Rstudio` or `Launch VSCodium`, but if you do that you will have access to _only_ that application instead of the Ubuntu desktop and all installed applications.

1. A new window should open in your browser with the Ubuntu desktop view.
You can now continue to [develop analysis using your virtual computer](STUB-LINK to developing with LSfR).
