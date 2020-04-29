{{- define "imagePullSecret" }}
{{- printf "{\"auths\": {\"%s\": {\"auth\": \"%s\"}}}" .Values.gitlab_registry_url (printf "%s:%s" .Values.gitlab_deploy_user .Values.gitlab_deploy_password | b64enc) | b64enc }}
{{- end }}