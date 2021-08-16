{{/* Generate a template that can be used for both assigned and unassigned xperiments */}}
{{- define "worker.pod-template" -}}
    metadata:
      labels:
        sandboxId: "{{ .Values.sandboxId }}"
    spec:
      containers:
      - name: "{{ .Release.Name }}-r"
        image: "{{ .Values.r.image }}"
        volumeMounts:
        - name: 'data'
          mountPath: '/data'
        - name: watch-script
          mountPath: /var/lib/watchfile
          readOnly: true
        - name: shutdown-file
          mountPath: /var/lib/shutdown-file
        - name: podinfo
          mountPath: /etc/podinfo
        ports:
        - containerPort: 4000
        resources:
          requests:
            memory: "27Gi"
      - name: "{{ .Release.Name }}"
        image: "{{ .Values.python.image }}"
        env:
        - name: AWS_XRAY_DAEMON_ADDRESS
          value: xray-service.default:2000
        - name: 'K8S_ENV'
          value: {{ .Values.kubernetes.env | quote }}
        - name: 'IGNORE_TIMEOUT'
          valueFrom:
            configMapKeyRef:
              name: instance-config
              key: ignoreTimeout
        volumeMounts:
        - name: 'data'
          mountPath: '/data'
        - name: watch-script
          mountPath: /var/lib/watchfile
          readOnly: true
        - name: shutdown-file
          mountPath: /var/lib/shutdown-file
        - name: podinfo
          mountPath: /etc/podinfo
        resources:
          requests:
            memory: "2Gi"
      volumes:
      - name: 'data'
      - name: watch-script
        configMap:
          name: "watch-script"
          items:
            - key: watcher.sh
              path: watcher.sh
            - key: entrypoint.sh
              path: entrypoint.sh
      - name: shutdown-file
      - name: podinfo
        downwardAPI:
          items:
            - path: "labels"
              fieldRef:
                fieldPath: metadata.labels
      restartPolicy: Always
      serviceAccountName: 'deployment-runner'
{{- end -}}