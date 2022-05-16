
VERSION := $(shell cat VERSION)
DOCKER_USERNAME := sbastkowski
DOCKER_IMAGENAME := albatradis
DOCKER_PATH := ${DOCKER_USERNAME}/${DOCKER_IMAGENAME}

dev:
	python3 setup.py develop

install:
	python3 setup.py install

unit_test:
	pytest --cov=albatradis --cov-report=xml --cov-branch --doctest-modules

script_test:
	albatradis -v -a data/albatradis_data/reference_BW25113_short.embl data/albatradis_data/025mgLTricRep1.insert_site_plot_short.gz  data/albatradis_data/controlLBrep1.insert_site_plot_short.gz
	albatradis-presence_absence data/presence_absence_data/reference_BW25113.embl data/presence_absence_data/gene_report_0008mgL.csv data/presence_absence_data/gene_report_0015mgL.csv data/presence_absence_data/gene_report_003mgL.csv data/presence_absence_data/gene_report_006mgL.csv data/presence_absence_data/gene_report_0125mgL.csv data/presence_absence_data/gene_report_025mgL.csv data/presence_absence_data/gene_report_05mgL.csv data/presence_absence_data/gene_report_1mgL.csv
	rm -r output

test: unit_test script_test

sonar: test
	sonar-scanner

clean:
	rm -rf build dist albatradis.egg-info coverage.xml .scannerwork

docker-build:
	docker build -t ${DOCKER_PATH}:${VERSION} .
	docker tag ${DOCKER_PATH}:${VERSION} ${DOCKER_PATH}:dev

docker-push: docker-build
	echo "${DOCKER_PASSWORD}" | docker login -u "${DOCKER_USERNAME}" --password-stdin
	docker tag ${DOCKER_PATH}:${VERSION} ${DOCKER_PATH}:latest
	docker push ${DOCKER_PATH}:latest
	docker push ${DOCKER_PATH}:${VERSION}


update_master_branch:
	git fetch origin
	git checkout master
	git pull origin master

bump_major_version: update_master_branch
	python3 bump_version.py --mode=major

bump_minor_version: update_master_branch
	python3 bump_version.py --mode=minor

bump_version: update_master_branch
	python3 bump_version.py --mode=patch

release:
	git commit VERSION -m "chore(package): Bump version up to $(shell cat VERSION)"
	git tag "$(shell cat VERSION)"
	git push origin master
	git push origin master --tags


release_major: bump_major_version release
	echo "Released new major version"

release_minor: bump_minor_version release
	echo "Released new minor version"

release_patch: bump_version release
	echo "Released new patch"


